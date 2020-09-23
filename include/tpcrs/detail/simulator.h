#pragma once

#include <vector>
#include <utility>

#include "TH1.h"
#include "TRandom.h"

#include "tpcrs/tpcrs_core.h"
#include "tpcrs/detail/coords.h"
#include "tpcrs/detail/distorter.h"
#include "tpcrs/detail/mag_field.h"
#include "tpcrs/detail/TF1F.h"
#include "tpcrs/detail/track_helix.h"


namespace tpcrs { namespace detail {


class Simulator
{
 public:

  Simulator(const tpcrs::Configurator& cfg);

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized) const;

  template<typename InputIt, typename OutputIt, typename MagField>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized, const MagField& mag_field) const;

  template<typename InputIt, typename OutputIt>
  OutputIt Distort(InputIt first_hit, InputIt last_hit, OutputIt distorted) const;

 private:

  template<typename InputIt, typename OutputIt1, typename OutputIt2, typename MagField>
  void Simulate(InputIt first_hit, InputIt last_hit, OutputIt1 digitized, OutputIt2 distorted, const MagField& mag_field) const;

  struct TrackSegment {
    int charge;
    double mass;
    tpcrs::SimulatedHit simu_hit;
    /// The original coordinates of the hit with applied distortions
    StTpcLocalSectorCoordinate coorLS;
    StTpcLocalSectorDirection  dirLS;
    StTpcLocalSectorDirection  BLS;
    /// Hardware coordinates corresponding to the local to the sector `coorLS`
    StTpcPadCoordinate Pad;
  };

  using ChargeContainer = std::vector<tpcrs::SimulatedCharge>;

  enum InOut {kInner = 0, kOuter = 1};

  enum class dEdxModel : unsigned {
    kBichsel = 1,
    kHeed = 2
  };

  enum {kPadMax = 32, kTimeBacketMax = 64};

  const tpcrs::Configurator& cfg_;
  const CoordTransform transform_;
  tpcrs::DigiChannelMap digi_;

  static double shapeEI(double* x, double* par = 0);
  static double shapeEI(double t, double t0, double tau_I, double tau_C);
  static double shapeEI_I(double* x, double* par = 0);
  static double shapeEI_I(double t, double timebin_width, double norm, int io);
  static double shapeEI3(double* x, double* par = 0);
  static double shapeEI3(double t, double t0, double tau_F, double tau_P, double tau_I, double tau_C);
  static double shapeEI3_I(double* x, double* par = 0);
  static double shapeEI3_I(double x, double timebin_width, double norm, int io);
  static double fei(double t, double t0, double T);
  static double polya(double* x, double* par);
  static double Ec(double* x, double* p); // minimal energy to create an ion pair
  static double PadResponseFunc(double* x, double* p);
  static double PadResponseFunc(double x, double w, double h, double K3, double cross_talk, double p);
  static double Gatti(double x, double pad_width, double anode_cathode_gap, double K3);
  static double InducedCharge(double s, double h, double ra, double Va, double &t0);
  static void ParticleProperties(int particle_id, int& charge, double& mass);

  /**
   * sigma = electron_range*(eEnery/electron_range_energy)^electron_range_power
   *
   * electron_range = 0.0055 cm
   * electron_range_energy = 3000 eV
   * electron_range_power = 1.78
   */
  static double ElectronRange(float dE, float dEr)
  {
    return dE > 3000. ? 0.0055 * std::pow((dE + dEr) / 3000., 1.78) : 0;
  }

  static double GatingGridTransparency(double x, double x_min = 10, double x_max = 56)
  {
    if (x <= x_min || x >= x_max) return 1;
    return std::max( 0., 1 - 6.27594134307865925e+00 * std::exp(-2.87987e-01 * (x - 1.46222e+01)) );
  };

  double GetNoPrimaryClusters(double betaGamma, int charge) const;

  template<typename OutputIt>
  void DigitizeSector(unsigned int sector, const ChargeContainer& binned_charge, OutputIt digitized) const;

  void SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail) const;
  void SimulateAsic(std::vector<short>& ADC) const;

  template<typename InputIt, typename OutputIt, typename MagField>
  TrackSegments CreateTrackSegments(InputIt first_hit, InputIt last_hit, OutputIt distorted, const MagField& mag_field) const;

  template<typename OutputIt, typename MagField>
  TrackSegment CreateTrackSegment(const tpcrs::SimulatedHit& hit, OutputIt distorted, const MagField& mag_field) const;

  double CalcBaseGain(int sector, int row) const;
  double CalcLocalGain(const TrackSegment& segment) const;

  void SignalFromSegment(const TrackSegment& segment, TrackHelix track,
    double gain_local,
    ChargeContainer& binned_charge, int& nP, double& dESum, double& dSSum) const;

  void LoopOverElectronsInCluster(
    const std::vector<float>& rs, const TrackSegment& segment, ChargeContainer& binned_charge,
    double xRange, Coords xyzC, double gain_local) const;

  void GenerateSignal(const TrackSegment &segment, Coords at_readout, int rowMin, int rowMax,
                      const TF1F* shaper, ChargeContainer& binned_charge, double gain_local_gas) const;

  std::vector<float> NumberOfElectronsInCluster(const TF1& heed, float dE, float& dEr) const;

  Coords TransportToReadout(const Coords c, double omega_tau, bool& missed_readout, bool& is_ground_wire) const;

  double dEdxCorrection(const TrackSegment &segment) const;

  using FuncParams_t = std::vector< std::pair<std::string, double> >;

  /// Initializes the mPadResponseFunction
  void InitPadResponseFuncs(int io, int sector);

  void InitChargeFractionFuncs(int io, int sector);

  /// Initializes the mShaperResponses array with shape functions
  void InitShaperFuncs(int io, int sector, std::array<std::vector<TF1F>, 2>& funcs,
    double (*shape)(double*, double*), FuncParams_t params, double timeBinMin, double timeBinMax);

  void InitAlphaGainVariations(double t0IO[2]);

  static TF1 fgTimeShape3[2];
  static TF1 fgTimeShape0[2];

  /// dEdx model
  dEdxModel dEdx_model_;

  /// The number of primary electrons per unit length as a function of
  /// beta*gamma
  TH1D  dNdx_;
  TH1D  dNdx_log10_;

  /// A probability distribution for energy lost by primary electrons.
  /// Can be parameterized as
  /// dE = TpcResponseSimulator.W * gRandom->Poisson(TpcResponseSimulator.Cluster);
  TH1D  dNdE_log10_;

  std::array<std::vector<TF1F>, 2>  mShaperResponses;
  std::vector<TF1F>  mChargeFraction;
  std::vector<TF1F>  mPadResponseFunction;
  std::vector<TF1F>  mPolya;
  TF1    mHeed;

  std::vector<double> alpha_gain_variations_;
};


template<typename InputIt, typename OutputIt>
OutputIt Simulator::Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized) const
{
  MagField mag_field(cfg_);
  return Digitize(first_hit, last_hit, digitized, mag_field);
}


template<typename InputIt, typename OutputIt, typename MagField>
OutputIt Simulator::Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized, const MagField& mag_field) const
{
  std::vector<tpcrs::DistortedHit> dummy;
  Simulate(first_hit, last_hit, digitized, std::back_inserter(dummy), mag_field);
  return digitized;
}


template<typename InputIt, typename OutputIt>
OutputIt Simulator::Distort(InputIt first_hit, InputIt last_hit, OutputIt distorted) const
{
  tpcrs::detail::MagField mag_field(cfg_);
  std::vector<tpcrs::DigiHit> dummy;
  Simulate(first_hit, last_hit, std::back_inserter(dummy), distorted, mag_field);
  return distorted;
}


template<typename InputIt, typename OutputIt1, typename OutputIt2, typename MagField>
void Simulator::Simulate(InputIt first_hit, InputIt last_hit, OutputIt1 digitized, OutputIt2 distorted, const MagField& mag_field) const
{
  std::vector<TrackSegment> segments = CreateTrackSegments(first_hit, last_hit, distorted, mag_field);

  static int nCalls = 0;
  gRandom->SetSeed(2345 + nCalls++);

  ChargeContainer binned_charge(digi_.total_timebins(), {0, 0});

  for (auto segment_iter = begin(segments); segment_iter != end(segments); ++segment_iter)
  {
    auto segment = *segment_iter;
    unsigned curr_sector =      segment_iter ->simu_hit.volume_id % 10000 / 100;
    unsigned next_sector = next(segment_iter)->simu_hit.volume_id % 10000 / 100;

    bool boundary = next_sector != curr_sector;

    if (segment.charge == 0 || segment.Pad.timeBucket < 0 || segment.Pad.timeBucket > digi_.n_timebins)
    {
      if (boundary) {
        DigitizeSector(curr_sector, binned_charge, digitized);
        ChargeContainer(digi_.total_timebins(), {0, 0}).swap(binned_charge);
      }
      continue;
    }

    // Calculate local gain corrected for dE/dx
    double gain_local = CalcLocalGain(segment);
    if (gain_local == 0)
    {
      if (boundary) {
        DigitizeSector(curr_sector, binned_charge, digitized);
        ChargeContainer(digi_.total_timebins(), {0, 0}).swap(binned_charge);
      }
      continue;
    }

    // Initialize propagation
    // Magnetic field BField must be in kilogauss
    // kilogauss = 1e-1*tesla = 1e-1*(volt*second/meter2) = 1e-1*(1e-6*1e-3*1/1e4) = 1e-14
    TrackHelix track(segment.dirLS.position,
                     segment.coorLS.position,
                     segment.BLS.position.z * 1e-14 * segment.charge, 1);
    // Propagate track to the middle of the pad row plane defined by the
    // nominal center point and the normal in this sector coordinate system
    double sR = track.pathLength({0, tpcrs::RadialDistanceAtRow(segment.Pad.row, cfg_), 0}, {0, 1, 0});

    // Update hit position based on the new track crossing the middle of pad row
    if (sR < 1e10) {
      segment.coorLS.position = {track.at(sR).x, track.at(sR).y, track.at(sR).z};
      transform_.local_sector_to_hardware(segment.coorLS, segment.Pad);
    }

    int nP = 0;
    double dESum = 0;
    double dSSum = 0;

    SignalFromSegment(segment, track, gain_local, binned_charge, nP, dESum, dSSum);

      *distorted = tpcrs::DistortedHit{
        {segment.coorLS.position.x, segment.coorLS.position.y, segment.coorLS.position.z},
        {segment.dirLS.position.x,  segment.dirLS.position.y,  segment.dirLS.position.z},
        dESum * 1e-9, // electronvolt in GeV
        dSSum,
        nP
      };

    if (boundary) {
      DigitizeSector(curr_sector, binned_charge, digitized);
      ChargeContainer(digi_.total_timebins(), {0, 0}).swap(binned_charge);
    }
  }
}


template<typename InputIt, typename OutputIt, typename MagField>
Simulator::TrackSegments Simulator::CreateTrackSegments(InputIt first_hit, InputIt last_hit, OutputIt distorted, const MagField& mag_field) const
{
  TrackSegments segments;
  for (auto curr_hit = first_hit; curr_hit != last_hit; ++curr_hit)
  {
    segments.push_back( CreateTrackSegment(*curr_hit, distorted, mag_field) );
  }
  return segments;
}


template<typename OutputIt, typename MagField>
Simulator::TrackSegment Simulator::CreateTrackSegment(const tpcrs::SimulatedHit& hit, OutputIt distorted, const MagField& mag_field) const
{
  int sector = hit.volume_id % 10000 / 100;

  TrackSegment segment{};
  StGlobalCoordinate xyzG{hit.x, hit.y, hit.z};
  segment.simu_hit = hit;
  ParticleProperties(hit.particle_id, segment.charge, segment.mass);

  StTpcLocalSectorCoordinate coorS;
  // GlobalCoord -> LocalSectorCoord. This transformation can result in a row
  // that is not the same as (volId % 100)
  transform_.global_to_local_sector(xyzG, coorS, sector, 0);
  StTpcLocalCoordinate coorLT;  // before distortions
  transform_.global_to_local(xyzG, coorLT, sector, coorS.row);

  // Get magnetic field at the hit position
  auto B_field = mag_field.ValueAt( Coords{hit.x, hit.y, hit.z});
  // distortion and misalignment
  // replace pxy => direction and try linear extrapolation
  Coords pxyzG{hit.px, hit.py, hit.pz};
  StGlobalDirection dirG{pxyzG.unit()}; // XXX Why scale the momentum?
  // TODO: Remove cast to float when new reference is introduced for tests
  StGlobalDirection BG{float(B_field.x), float(B_field.y), float(B_field.z)};
  transform_.global_to_local_sector_dir( dirG, segment.dirLS, sector, coorS.row);
  transform_.global_to_local_sector_dir(   BG, segment.BLS,   sector, coorS.row);

  // Distortions
  static Distorter distorter(cfg_);
  coorLT.position = distorter.Distort(coorLT.position, coorLT.sector, mag_field);
  transform_.local_to_global(coorLT, xyzG);

  transform_.local_to_local_sector(coorLT, segment.coorLS);

  *distorted = tpcrs::DistortedHit{
    xyzG.position.x, xyzG.position.y, xyzG.position.z,
    dirG.position.x, dirG.position.y, dirG.position.z
  };

  double driftLength = segment.coorLS.position.z + hit.tof * tpcrs::DriftVelocity(sector, cfg_);

  if (driftLength > -1.0 && driftLength <= 0) {
    if ((!tpcrs::IsInner(coorS.row, cfg_) && driftLength > - cfg_.S<tpcWirePlanes>().outerSectorAnodeWirePadSep) ||
        ( tpcrs::IsInner(coorS.row, cfg_) && driftLength > - cfg_.S<tpcWirePlanes>().innerSectorAnodeWirePadSep))
      driftLength = std::abs(driftLength);
  }

  segment.coorLS.position.z = driftLength;
  transform_.local_sector_to_hardware(segment.coorLS, segment.Pad);

  return segment;
}


template<typename OutputIt>
void Simulator::DigitizeSector(unsigned int sector, const ChargeContainer& binned_charge, OutputIt digitized) const
{
  double pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX;
  double ped = cfg_.S<TpcResponseSimulator>().AveragePedestal;

  std::vector<short> ADCs_(binned_charge.size(), 0);

  auto bc = binned_charge.begin();
  auto adcs_iter = ADCs_.begin();

  for (auto ch = digi_.channels.begin(); ch != digi_.channels.end(); ch += digi_.n_timebins)
  {
    double gain = cfg_.S<tpcPadGainT0>().Gain[sector-1][ch->row-1][ch->pad-1];

    if (gain <= 0) {
      bc        += digi_.n_timebins;
      adcs_iter += digi_.n_timebins;
      continue;
    }

    for (int i=0; i != digi_.n_timebins; ++i, ++bc, ++adcs_iter)
    {
      int adc = int(bc->Sum / gain + gRandom->Gaus(ped, pedRMS) - ped);
      // Zero negative values
      adc = adc & ~(adc >> 31);
      // Select minimum between adc and 1023, i.e. overflow at 1023
      adc = adc - !(((adc - 1023) >> 31) & 0x1) * (adc - 1023);

      *adcs_iter = adc;
    }
  }

  for (auto adcs_iter = ADCs_.begin(); adcs_iter != ADCs_.end(); adcs_iter += digi_.n_timebins)
  {
    SimulateAltro(adcs_iter, adcs_iter + digi_.n_timebins, true);
  }

  auto ch = digi_.channels.begin();
  adcs_iter = ADCs_.begin();

  for (auto bc = binned_charge.begin(); bc != binned_charge.end(); ++bc, ++ch, ++adcs_iter)
  {
    if (*adcs_iter == 0) continue;
    *digitized = tpcrs::DigiHit{sector, ch->row, ch->pad, ch->timebin, *adcs_iter, bc->TrackId};
  }
}

} }
