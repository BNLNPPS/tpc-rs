#ifndef TPCRS_SIMULATOR_H_
#define TPCRS_SIMULATOR_H_

#include <vector>
#include <utility>

#include "TF1F.h"
#include "TH1.h"

#include "tpcrs/structs.h"
#include "tpcrs/tpcrs.h"
#include "coords.h"
#include "dedx_correction.h"
#include "struct_containers.h"
#include "track_helix.h"
#include "mag_field_utils.h"


class Simulator
{
 public:

  Simulator(const tpcrs::Configurator& cfg);
  ~Simulator();

  template<typename InputIt, typename OutputIt1, typename OutputIt2>
  void Simulate(InputIt first_hit, InputIt last_hit, OutputIt1 digi_data, OutputIt2);

 private:

  using SimuHitIt = std::vector<tpcrs::SimulatedHit>::iterator;
  using DigiInserter = std::back_insert_iterator<std::vector<tpcrs::DigiHit>>;
  using DistInserter = std::back_insert_iterator<std::vector<tpcrs::DistortedHit>>;

  struct TrackSegment {
    int charge;
    double mass;
    tpcrs::SimulatedHit* simu_hit;
    /// The original coordinates of the hit with applied distortions
    StTpcLocalSectorCoordinate coorLS;
    StTpcLocalSectorDirection  dirLS;
    StTpcLocalSectorDirection  BLS;
    /// Hardware coordinates corresponding to the local to the sector `coorLS`
    StTpcPadCoordinate Pad;
  };

  using ChargeContainer = std::vector<tpcrs::SimulatedCharge>;

  enum InOut {kInner = 0, kOuter = 1};
  enum EMode {kPAI         = 0,// switch to PAI from GEANT (obsolete)
              kBICHSEL     = 1,// switch to Bichsel from GEANT
              kHEED        = 6,// switch to HEED
              kGAINOAtALL  = 2,// do not use GAIN at all
              kdEdxCorr    = 3,// do use TpcdEdxCorrection
              kDistortion  = 4,// include distortions
              kNoToflight  = 5 // don't account for particle time of flight
             };
  enum {kPadMax = 32, kTimeBacketMax = 64, kRowMax = 72};

  const tpcrs::Configurator& cfg_;
  const CoordTransform transform_;
  MagFieldUtils mag_field_utils_;
  const int num_sectors_;            //!
  const int max_rows_;               //!
  const int max_pads_;               //!
  const int max_timebins_;

  static double shapeEI(double* x, double* par = 0);
  static double shapeEI_I(double* x, double* par = 0);
  static double shapeEI3(double* x, double* par = 0);
  static double shapeEI3_I(double* x, double* par = 0);
  static double fei(double t, double t0, double T);
  static double polya(double* x, double* par);
  static double Ec(double* x, double* p); // minimal energy to create an ion pair
  static double PadResponseFunc(double* x, double* p);
  static double Gatti(double* x, double* p);
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

  bool IsInner(int row, int sector=1)
  {
    return row <= cfg_.C<St_tpcPadConfigC>().numberOfInnerRows(sector);
  };

  double GetNoPrimaryClusters(double betaGamma, int charge);
  void DigitizeSector(unsigned int sector, const ChargeContainer& binned_charge, DigiInserter digi_data);

  /// Returns the number of non-zero elements in ADC
  void SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail);
  void SimulateAsic(std::vector<short>& ADC);

  void CreateTrackSegments(SimuHitIt, SimuHitIt, std::vector<TrackSegment>&);

  TrackSegment CreateTrackSegment(tpcrs::SimulatedHit& hit);

  double CalcBaseGain(int sector, int row);
  double CalcLocalGain(int sector, int row, double gain_base, double dedx_corr);

  void SignalFromSegment(const TrackSegment& segment, TrackHelix track,
    double gain_local,
    ChargeContainer& binned_charge, int& nP, double& dESum, double& dSSum);

  void LoopOverElectronsInCluster(
    std::vector<float> rs, const TrackSegment& segment, ChargeContainer& binned_charge,
    double xRange, Coords xyzC, double gain_local);

  void GenerateSignal(const TrackSegment &segment, int rowMin, int rowMax,
                      TF1F* shaper, ChargeContainer& binned_charge, double gain_local_gas);

  std::vector<float> NumberOfElectronsInCluster(const TF1& heed, float dE, float& dEr);

  double dEdxCorrection(const TrackSegment &segment);

  using FuncParams_t = std::vector< std::pair<std::string, double> >;

  /// Initializes the mShaperResponses array with shape functions
  void InitShaperFuncs(int io, int sector, std::array<std::vector<TF1F>, 2>& funcs,
    double (*shape)(double*, double*), FuncParams_t params, double timeBinMin, double timeBinMax);

  static TF1F     fgTimeShape3[2];  //!
  static TF1F     fgTimeShape0[2];   //!
  int    options_;
  TH1D*    mdNdx;                     //!
  TH1D*    mdNdxL10;                  //!
  TH1D*    mdNdEL10;                  //!
  std::array<std::vector<TF1F>, 2>  mShaperResponses;     //!
  std::array<std::vector<TF1F>, 2>  mChargeFraction;      //!
  std::array<std::vector<TF1F>, 2>  mPadResponseFunction; //!
  TF1F   mPolya[2];                   //!
  TF1    mHeed;                       //!
  StTpcdEdxCorrection dEdx_correction_; // !
  double InnerAlphaVariation[24];   //!
  double OuterAlphaVariation[24];   //!
  int    numberOfInnerSectorAnodeWires; //!
  double firstInnerSectorAnodeWire; //!
  double lastInnerSectorAnodeWire;  //!
  int    numberOfOuterSectorAnodeWires; //!
  double firstOuterSectorAnodeWire; //!
  double lastOuterSectorAnodeWire;  //!
  double anodeWirePitch;            //!
  double anodeWireRadius;           //!
  double innerSectorAnodeVoltage[24];//!
  double outerSectorAnodeVoltage[24];//!
  double xOnWire, yOnWire, zOnWire; //!

  tpcrs::DigiChannelMap digi_;
};


template<>
void Simulator::Simulate(SimuHitIt first_hit, SimuHitIt last_hit, DigiInserter digi_data, DistInserter);

#endif
