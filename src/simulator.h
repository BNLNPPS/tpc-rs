#ifndef TPCRS_SIMULATOR_H_
#define TPCRS_SIMULATOR_H_

#include <vector>
#include <utility>

#include "TF1F.h"
#include "TH1.h"

#include "tpcrs/digi_data.h"
#include "tpcrs/structs.h"
#include "tpc_db.h"

class Altro;
class StTpcdEdxCorrection;
class HitPoint_t;
class SignalSum_t;

class StTpcRSMaker
{
 public:

  StTpcRSMaker(double eCutOff = 1e-3, const char* name = "TpcRS");
  ~StTpcRSMaker();

  void Make(const std::vector<g2t_tpc_hit_st>& g2t_tpc_hit,
            const std::vector<g2t_track_st>& g2t_track,
            const std::vector<g2t_vertex_st>& g2t_vertex, tpcrs::DigiData& digi_data);

 private:

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
  void DigitizeSector(int sector, tpcrs::DigiData& digi_data, std::vector<SignalSum_t>& SignalSum);
  int    AsicThresholds(short* ADCs);
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

  void BuildTrackSegments(int sector, const std::vector<size_t>& sorted_index, int sortedIndex,
    const g2t_tpc_hit_st* tpc_hit_begin, const g2t_vertex_st& geant_vertex,
    HitPoint_t TrackSegmentHits[100], double& smin, double& smax, int& nSegHits, int& sIndex);

  void TrackSegment2Propagate(g2t_tpc_hit_st* tpc_hitC, const g2t_vertex_st& geant_vertex, HitPoint_t &TrackSegmentHits, double& smin, double& smax);
  void GenerateSignal(const HitPoint_t &TrackSegmentHits, int sector, int rowMin, int rowMax, double sigmaJitterT,
                      double sigmaJitterX, TF1F* shaper, std::vector<SignalSum_t>& SignalSum,
                      double& total_signal_in_cluster, double gain_local, double gain_gas);

  double dEdxCorrection(HitPoint_t &TrackSegmentHits);

  static double GatingGridTransparency(double x, double x_min = 10, double x_max = 56)
  {
    if (x <= x_min || x >= x_max) return 1;
    return std::max( 0., 1 - 6.27594134307865925e+00 * std::exp(-2.87987e-01 * (x - 1.46222e+01)) );
  };

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
  std::array<std::array<TF1F*, 24>, 2>  mChargeFraction;      //!
  std::array<std::array<TF1F*, 24>, 2>  mPadResponseFunction; //!
  TF1F*  mPolya[2];                   //!
  TF1*   mHeed;                       //!
  StTpcdEdxCorrection* m_TpcdEdxCorrection; // !
  double InnerAlphaVariation[24];   //!
  double OuterAlphaVariation[24];   //!
  Altro* mAltro;                      //!
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
  double    mLocalYDirectionCoupling[2][24][7]; //!
  double xOnWire, yOnWire, zOnWire; //!
  int pad0;                         //!
  int tbk0;                         //!
  double padsdE[kPadMax];           //!
  double tbksdE[kTimeBacketMax];    //!
  double rowsdEH[kRowMax];          //!
  double rowsdE[kRowMax];           //!
  const double min_signal_;           //!
  const double electron_range_;       //!
  const double electron_range_energy_; //!
  const double electron_range_power_;  //!
  const double max_electron_energy_;                 //! cut for delta electrons
  const int max_sectors_;            //!
  const int max_pads_;               //!
  const int max_timebins_;
};
#endif
