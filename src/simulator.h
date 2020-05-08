#ifndef TPCRS_SIMULATOR_H_
#define TPCRS_SIMULATOR_H_

#include <vector>
#include <utility>

#include "TF1F.h"
#include "TH1.h"

#include "tpcrs/digi_data.h"
#include "tpcrs/structs.h"
#include "coords.h"
#include "dedx_correction.h"
#include "tpc_db.h"

class HitPoint_t;
class SignalSum_t;

class StTpcRSMaker
{
 public:

  StTpcRSMaker(double eCutOff = 1e-3, const char* name = "TpcRS");
  ~StTpcRSMaker();

  void Make(const std::vector<g2t_tpc_hit>& geant_hits,
            const std::vector<g2t_track>& geant_particles,
            const std::vector<g2t_vertex>& geant_vertices, tpcrs::DigiData& digi_data);

 private:

  const double min_signal_;           //!
  const double electron_range_;       //!
  const double electron_range_energy_; //!
  const double electron_range_power_;  //!
  const double max_electron_energy_;                 //! cut for delta electrons
  const int max_sectors_;            //!
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

  static double GatingGridTransparency(double x, double x_min = 10, double x_max = 56)
  {
    if (x <= x_min || x >= x_max) return 1;
    return std::max( 0., 1 - 6.27594134307865925e+00 * std::exp(-2.87987e-01 * (x - 1.46222e+01)) );
  };

  enum InOut {kInner, kOuter};
  enum EMode {kPAI         = 0,// switch to PAI from GEANT (obsolete)
              kBICHSEL     = 1,// switch to Bichsel from GEANT
              kHEED        = 6,// switch to HEED
              kGAINOAtALL  = 2,// do not use GAIN at all
              kdEdxCorr    = 3,// do use TpcdEdxCorrection
              kDistortion  = 4,// include distortions
              kNoToflight  = 5 // don't account for particle time of flight
             };
  enum {kPadMax = 32, kTimeBacketMax = 64, kRowMax = 72};
  int         Debug() const {return 1;}
  double GetNoPrimaryClusters(double betaGamma, int charge);
  void Print(Option_t* option = "") const;
  void DigitizeSector(int sector, tpcrs::DigiData& digi_data, std::vector<SignalSum_t>& binned_charge);
  int    AsicThresholds(short* ADCs);

  void BuildTrackSegments(int sector, const std::vector<size_t>& sorted_index, int sortedIndex,
    const g2t_tpc_hit* tpc_hit_begin, const g2t_vertex& geant_vertex,
    HitPoint_t TrackSegmentHits[100], double& smin, double& smax, int& nSegHits, int& sIndex);

  void TrackSegment2Propagate(g2t_tpc_hit* tpc_hitC, const g2t_vertex& geant_vertex, HitPoint_t &TrackSegmentHits, double& smin, double& smax);

  double LoopOverElectronsInCluster(std::vector<float> rs, const HitPoint_t& TrackSegmentHits, std::vector<SignalSum_t>& binned_charge,
    int sector, int row,
    double xRange, Coords xyzC, double gain_local,
    double SigmaT, double SigmaL, double OmegaTau);

  void GenerateSignal(const HitPoint_t &TrackSegmentHits, int sector, int rowMin, int rowMax, double sigmaJitterT,
                      TF1F* shaper, std::vector<SignalSum_t>& binned_charge,
                      double& total_signal_in_cluster, double gain_local, double gain_gas);

  std::vector<float> NumberOfElectronsInCluster(const TF1& heed, float dE, float& dEr);

  double dEdxCorrection(const HitPoint_t &path_segment);

  using FuncParams_t = std::vector< std::pair<std::string, double> >;

  /// Initializes the mShaperResponses array with shape functions
  void InitShaperFuncs(int io, int sector, std::array<std::array<TF1F*, 24>, 2>& funcs,
    double (*shape)(double*, double*), FuncParams_t params, double timeBinMin, double timeBinMax);

  static TF1F*     fgTimeShape3[2];  //!
  static TF1F*     fgTimeShape0[2];   //!
  int    options_;
  TH1D*    mdNdx;                     //!
  TH1D*    mdNdxL10;                  //!
  TH1D*    mdNdEL10;                  //!
  std::array<std::array<TF1F*, 24>, 2>  mShaperResponses;     //!
  std::array<std::vector<TF1F>, 2>  mChargeFraction;      //!
  std::array<std::vector<TF1F>, 2>  mPadResponseFunction; //!
  TF1F   mPolya[2];                   //!
  TF1    mHeed;                       //!
  StTpcdEdxCorrection m_TpcdEdxCorrection; // !
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
  double rowsdE[kRowMax];           //!
};
#endif
