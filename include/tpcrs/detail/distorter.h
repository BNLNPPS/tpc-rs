#pragma once

#include "tpcrs/configurator.h"
#include "tpcrs/detail/enum_bitset.h"
#include "struct_containers.h"


namespace tpcrs { namespace detail {


struct DistorterConfigurator
{
  DistorterConfigurator(const tpcrs::Configurator& cfg, double mag_field_z = 0)
  {
    double cathode_voltage    = cfg.S<tpcHighVoltages>().cathode * 1000;
    double gated_grid_voltage = cfg.S<tpcHighVoltages>().gatedGridRef;

    double tensorV1 = cfg.S<tpcOmegaTau>().tensorV1;
    double tensorV2 = cfg.S<tpcOmegaTau>().tensorV2;

    double drift_velocity = 1e-6 * tpcrs::DriftVelocity(24, cfg);

    readout_plane_z = cfg.S<tpcPadPlanes>().outerSectorPadPlaneZ - cfg.S<tpcWirePlanes>().outerSectorGatingGridPadSep;

    // Theoretically, omega_tau is defined as
    //
    // omega_tau = -10. * B.z * drift_velocity / electric_field;  // cm/microsec, Volts/cm
    //
    // We use definitions from Amendolia et al NIM A235 (1986) 296 and include
    // their characterization of the electron drift velocity tensor with
    // different omega-tau's in different directions.
    //
    // tensorV1 is the drift velocity tensor term in the ExB direction
    // tensorV2 is the drift velocity tensor term in the direction perpendicular
    // to Z and ExB

    // Electric Field (V/cm) Magnitude
    double electric_field = std::abs((cathode_voltage - gated_grid_voltage) / readout_plane_z);

    // For an electron, omega_tau carries the sign opposite of B (mag_field_z)
    // B is in kGauss, the sign of B is important
    omega_tau = -10.0 * mag_field_z * drift_velocity / electric_field;

    const_0 = 1. / (1. +  tensorV2 * tensorV2 * omega_tau * omega_tau);
    const_1 = tensorV1 * omega_tau / (1. + tensorV1 * tensorV1 * omega_tau * omega_tau);
    const_2 = tensorV2 * tensorV2 * omega_tau * omega_tau / (1. + tensorV2 * tensorV2 * omega_tau * omega_tau);
  }

  // Z location of TPC readout plane (absolute value in cm)
  double readout_plane_z;
  double omega_tau;
  double const_0;
  double const_1;
  double const_2;
};


class Distorter
{
 public:

  enum class Distortions : unsigned {
    kNone = 0,
    kMagneticField
  };

  Distorter(const tpcrs::Configurator& cfg, double mag_field_z, Distortions distortions = Distortions::kNone) :
    cfg_(cfg),
    distortions_(distortions),
    dcfg_(cfg, mag_field_z)
  {}

  template<typename Vec3, typename MagField>
  Vec3 Distort(Vec3 p, int sector, const MagField& mag_field) const;

 private:

  const tpcrs::Configurator& cfg_;

  BitSet<Distortions> distortions_;

  DistorterConfigurator dcfg_;

  /** Distortion due to the shape of the magnetic field */
  template<typename Vec3, typename MagField>
  Vec3 DistortDueMagneticField(Vec3 p, int sector, const MagField& mag_field) const;
};


template<typename Vec3, typename MagField>
Vec3 Distorter::Distort(Vec3 p, int sector, const MagField& mag_field) const
{
  Vec3 p_prime = p;

  if ( distortions_.test(Distortions::kMagneticField) ) {
    p_prime = DistortDueMagneticField(p, sector, mag_field);
  }

  return p_prime;
}


/**
 * Integrate backwards from TPC plane to the point the electron cluster was
 * born.
 */
template<typename Vec3, typename MagField>
Vec3 Distorter::DistortDueMagneticField(Vec3 p, int sector, const MagField& mag_field) const
{
  int sign = DetectorSide(sector, cfg_) == TPC::Half::first ? +1 : -1;

  // Prepare for different readout planes
  Vec3 p_prime{p.x, p.y, sign * dcfg_.readout_plane_z};

  // Determine the number of steps and ah (to be about 1.0 cm), n_steps must be
  // odd. ah carries the sign opposite of E (for forward integration)
  double ah;
  int n_steps;
  for (n_steps = 5 ; n_steps < 1000 ; n_steps += 2 ) {
    // Going Backwards! See note above.
    ah = (p.z - sign * dcfg_.readout_plane_z) / (n_steps - 1);

    if (std::abs(ah) < 1) break;
  }

  // Simpson's integration loop
  int flag = 1;

  for (int i = 1; i <= n_steps; ++i) {
    if (i == n_steps) flag = 1 ;

    p_prime.z += flag * (ah/3);

    // Work in kGauss, cm, use Cartesian coordinates
    auto B = mag_field.ValueAt(p_prime);

    if (std::abs(B.z) > 0.001) {
      p_prime.x += flag * (ah/3) * (dcfg_.const_2 * B.x - dcfg_.const_1 * B.y) / B.z;
      p_prime.y += flag * (ah/3) * (dcfg_.const_2 * B.y + dcfg_.const_1 * B.x) / B.z;
    }

    flag = (flag != 4 ? 4 : 2);
  }

  return p_prime;
}

} }
