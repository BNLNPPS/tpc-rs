#pragma once

#include <vector>
#include <utility>

#include "TF1F.h"
#include "TH1.h"
#include "TRandom.h"

#include "tpcrs/tpcrs_core.h"
#include "coords.h"
#include "dedx_correction.h"
#include "struct_containers.h"
#include "track_helix.h"
#include "mag_field_utils.h"


namespace tpcrs { namespace detail {


class Simulator
{
 public:

  Simulator(const tpcrs::Configurator& cfg);

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized);

  template<typename InputIt, typename OutputIt, typename MagField>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized, const MagField& mag_field);

  template<typename InputIt, typename OutputIt>
  OutputIt Distort(InputIt first_hit, InputIt last_hit, OutputIt distorted);

 private:

  template<typename InputIt, typename OutputIt1, typename OutputIt2>
  void Simulate(InputIt first_hit, InputIt last_hit, OutputIt1 digitized, OutputIt2 distorted);

  struct TrackSegment {
    int charge;
    double mass;
    const tpcrs::SimulatedHit* simu_hit;
    /// The original coordinates of the hit with applied distortions
    StTpcLocalSectorCoordinate coorLS;
    StTpcLocalSectorDirection  dirLS;
    StTpcLocalSectorDirection  BLS;
    /// Hardware coordinates corresponding to the local to the sector `coorLS`
    StTpcPadCoordinate Pad;
  };

  using ChargeContainer = std::vector<tpcrs::SimulatedCharge>;

  enum InOut {kInner = 0, kOuter = 1};
  enum EMode {kBICHSEL     = 1, /// Use Bichsel dE/dx model
              kHEED        = 6, /// Use Heed dE/dx model
              kGAINOAtALL  = 2, /// Do not use GAIN at all
              kdEdxCorr    = 3, /// Do use TpcdEdxCorrection
              kDistortion  = 4, /// Apply distortions
              kNoToflight  = 5  /// Do not account for particle time of flight
             };
  enum {kPadMax = 32, kTimeBacketMax = 64};

  const tpcrs::Configurator& cfg_;
  const CoordTransform transform_;
  tpcrs::DigiChannelMap digi_;
  MagFieldUtils mag_field_utils_;

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

  double GetNoPrimaryClusters(double betaGamma, int charge);

  template<typename OutputIt>
  void DigitizeSector(unsigned int sector, const ChargeContainer& binned_charge, OutputIt digitized);

  void SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail);
  void SimulateAsic(std::vector<short>& ADC);

  template<typename InputIt>
  void CreateTrackSegments(InputIt, InputIt, std::vector<TrackSegment>&);

  TrackSegment CreateTrackSegment(tpcrs::SimulatedHit& hit);

  double CalcBaseGain(int sector, int row);
  double CalcLocalGain(const TrackSegment& segment);

  void SignalFromSegment(const TrackSegment& segment, TrackHelix track,
    double gain_local,
    ChargeContainer& binned_charge, int& nP, double& dESum, double& dSSum);

  void LoopOverElectronsInCluster(
    std::vector<float> rs, const TrackSegment& segment, ChargeContainer& binned_charge,
    double xRange, Coords xyzC, double gain_local);

  void GenerateSignal(const TrackSegment &segment, Coords at_readout, int rowMin, int rowMax,
                      TF1F* shaper, ChargeContainer& binned_charge, double gain_local_gas);

  std::vector<float> NumberOfElectronsInCluster(const TF1& heed, float dE, float& dEr);

  Coords TransportToReadout(const Coords c, double OmegaTau, bool& missed_readout, bool& is_ground_wire);

  double dEdxCorrection(const TrackSegment &segment);

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
  int    options_;

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
  StTpcdEdxCorrection dEdx_correction_;

  std::vector<double> alpha_gain_variations_;
};


template<typename InputIt, typename OutputIt>
OutputIt Simulator::Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized)
{
  MagField mag_field(cfg_);
  return Digitize(first_hit, last_hit, digitized, mag_field);
}


template<typename InputIt, typename OutputIt, typename MagField>
OutputIt Simulator::Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized, const MagField& mag_field)
{
  std::vector<tpcrs::DistortedHit> dummy;
  Simulate(first_hit, last_hit, digitized, std::back_inserter(dummy));
  return digitized;
}


template<typename InputIt, typename OutputIt>
OutputIt Simulator::Distort(InputIt first_hit, InputIt last_hit, OutputIt distorted)
{
  std::vector<tpcrs::DigiHit> dummy;
  Simulate(first_hit, last_hit, std::back_inserter(dummy), distorted);
  return distorted;
}


template<typename InputIt, typename OutputIt1, typename OutputIt2>
void Simulator::Simulate(InputIt first_hit, InputIt last_hit, OutputIt1 digitized, OutputIt2 distorted)
{
  static int nCalls = 0;
  gRandom->SetSeed(2345 + nCalls++);

  std::vector< std::vector<TrackSegment> > segments_by_sector(digi_.n_sectors);
  std::vector<TrackSegment> segments_in_sector;

  auto first_hit_on_track = first_hit;
  int curr_direction = 0; // 0 - increase no of row, 1 - decrease no of. row.

  for (auto curr_hit = first_hit; curr_hit != last_hit; ++curr_hit)
  {
    auto next_hit = next(curr_hit);

    int  next_direction  = next_hit->volume_id % 100 - curr_hit->volume_id % 100 >= 0 ? 0 : 1;
    bool sector_boundary = next_hit->volume_id % 10000 / 100 !=
                           curr_hit->volume_id % 10000 / 100;
    bool track_boundary  = (next_hit->track_id != curr_hit->track_id || curr_hit->track_id == 0);

    bool start_new_track = track_boundary || curr_direction != next_direction || sector_boundary;

    if (start_new_track || next_hit == last_hit) {
      CreateTrackSegments(first_hit_on_track, next_hit, segments_in_sector);
      first_hit_on_track = next_hit;

      if ( (sector_boundary || next_hit == last_hit ) && segments_in_sector.size() != 0) {
        segments_by_sector[(curr_hit->volume_id % 10000) / 100 - 1] = segments_in_sector;
        segments_in_sector.clear();
      }
    }

    curr_direction = next_direction;
  }

  unsigned int sector = 1;
  for (auto segments_in_sector : segments_by_sector) {
    int nHitsInTheSector = 0;
    ChargeContainer binned_charge(digi_.total_timebins(), {0, 0});

    for (TrackSegment& segment : segments_in_sector) {

      // Calculate local gain corrected for dE/dx
      double gain_local = CalcLocalGain(segment);
      if (gain_local == 0) continue;

      // Initialize propagation
      // Magnetic field BField must be in kilogauss
      // kilogauss = 1e-1*tesla = 1e-1*(volt*second/meter2) = 1e-1*(1e-6*1e-3*1/1e4) = 1e-14
      TrackHelix track(segment.dirLS.position,
                       segment.coorLS.position,
                       segment.BLS.position.z * 1e-14 * segment.charge, 1);
      // Propagate track to middle of the pad row plane by the nominal center point and the normal
      // in this sector coordinate system
      double sR = track.pathLength({0, transform_.yFromRow(segment.Pad.sector, segment.Pad.row), 0}, {0, 1, 0});

      // Update hit position based on the new track crossing the middle of pad row
      if (sR < 1e10) {
        segment.coorLS.position = {track.at(sR).x, track.at(sR).y, track.at(sR).z};
        transform_.local_sector_to_hardware(segment.coorLS, segment.Pad, false, false); // don't use T0, don't use Tau
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

      nHitsInTheSector++;
    } // end do loop over segments for a given particle

    if (nHitsInTheSector) {
      DigitizeSector(sector, binned_charge, digitized);
    }

    sector++;
  }
}


template<typename InputIt>
void Simulator::CreateTrackSegments(InputIt first_hit, InputIt last_hit, std::vector<TrackSegment>& segments)
{
  for (auto ihit = first_hit; ihit != last_hit; ++ihit)
  {
    TrackSegment curr_segment = CreateTrackSegment(*ihit);

    if (curr_segment.charge == 0) continue;
    if (curr_segment.Pad.timeBucket < 0 || curr_segment.Pad.timeBucket > digi_.n_timebins) continue;

    segments.push_back(curr_segment);
  }
}


Simulator::TrackSegment Simulator::CreateTrackSegment(tpcrs::SimulatedHit& hit)
{
  int sector = hit.volume_id % 10000 / 100;

  TrackSegment segment;
  StGlobalCoordinate xyzG{hit.x[0], hit.x[1], hit.x[2]};
  segment.simu_hit = &hit;
  ParticleProperties(hit.particle_id, segment.charge, segment.mass);

  StTpcLocalSectorCoordinate coorS;
  // GlobalCoord -> LocalSectorCoord. This transformation can result in a row
  // that is not the same as (volId % 100)
  transform_.global_to_local_sector(xyzG, coorS, sector, 0);
  StTpcLocalCoordinate coorLT;  // before distortions
  transform_.global_to_local(xyzG, coorLT, sector, coorS.row);

  // Get magnetic field at the hit position
  float B_field[3];
  mag_field_utils_.GetFieldValue(hit.x, B_field);
  // distortion and misalignment
  // replace pxy => direction and try linear extrapolation
  Coords pxyzG{hit.p[0], hit.p[1], hit.p[2]};
  StGlobalDirection dirG{pxyzG.unit()};
  StGlobalDirection BG{B_field[0], B_field[1], B_field[2]};
  transform_.global_to_local_sector_dir( dirG, segment.dirLS, sector, coorS.row);
  transform_.global_to_local_sector_dir(   BG, segment.BLS,   sector, coorS.row);

  // Distortions
  if (TESTBIT(options_, kDistortion)) {
    float pos[3] = {(float ) coorLT.position.x, (float ) coorLT.position.y, (float ) coorLT.position.z};
    float posMoved[3];
    mag_field_utils_.DoDistortion(pos, posMoved, sector); // input pos[], returns posMoved[]
    coorLT.position = {posMoved[0], posMoved[1], posMoved[2]};       // after distortions
    transform_.local_to_global(coorLT, xyzG);
  }

  transform_.local_to_local_sector(coorLT, segment.coorLS);

  double driftLength = segment.coorLS.position.z + hit.tof * tpcrs::DriftVelocity(sector, cfg_);

  if (driftLength > -1.0 && driftLength <= 0) {
    if ((!tpcrs::IsInner(coorS.row, cfg_) && driftLength > - cfg_.S<tpcWirePlanes>().outerSectorAnodeWirePadSep) ||
        ( tpcrs::IsInner(coorS.row, cfg_) && driftLength > - cfg_.S<tpcWirePlanes>().innerSectorAnodeWirePadSep))
      driftLength = std::abs(driftLength);
  }

  segment.coorLS.position.z = driftLength;
  transform_.local_sector_to_hardware(segment.coorLS, segment.Pad, false, false); // don't use T0, don't use Tau

  return segment;
}


template<typename OutputIt>
void Simulator::DigitizeSector(unsigned int sector, const ChargeContainer& binned_charge, OutputIt digitized)
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
