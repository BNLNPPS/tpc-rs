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

  template<typename Vec3>
  Vec3 Distort(Vec3 p) const;

 private:

  const tpcrs::Configurator& cfg_;

  BitSet<Distortions> distortions_;

  DistorterConfigurator dcfg_;

  /** Distortion due to the shape of the magnetic field */
  template<typename Vec3>
  Vec3 DistortDueMagneticField(Vec3 p) const;
};



template<typename Vec3>
Vec3 Distorter::Distort(Vec3 p) const
{
  Vec3 p_prime{};

  if ( distortions_.test(Distortions::kMagneticField) ) {
    p_prime = DistortDueMagneticField(p);
  }

  return p_prime;
}


template<typename Vec3>
Vec3 Distorter::DistortDueMagneticField(Vec3 p) const
{
  return p;
}

} }
