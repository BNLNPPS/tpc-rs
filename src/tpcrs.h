#ifndef STAR_ST_TRS_MAKER_HH
#define STAR_ST_TRS_MAKER_HH

#include "StTpcRSMaker/TF1F.h"
#include "TH1.h"
#include "StTpcDb/StTpcDb.h"
#include "tpcrs/digi_data.h"

class Altro;
class StTpcdEdxCorrection;
class StTpcDigitalSector;
class HitPoint_t;
class g2t_tpc_hit_st;
class g2t_vertex_st;
class StTpcCoordinateTransform;

class St_g2t_tpc_hit;
class St_g2t_track;
class St_g2t_vertex;

struct SignalSum_t {
  float      Sum;
  short      Adc;
  short  TrackId;
};

class StTpcRSMaker
{
 public:
  enum EMode {kPAI         = 0,// switch to PAI from GEANT (obsolete)
              kBICHSEL     = 1,// switch to Bichsel from GEANT
              kHEED        = 6,// switch to HEED
              kGAINOAtALL  = 2,// do not use GAIN at all
              kdEdxCorr    = 3,// do use TpcdEdxCorrection
              kDistortion  = 4,// include distortions
              kNoToflight  = 5 // don't account for particle time of flight
             };
  enum {kPadMax = 32, kTimeBacketMax = 64, kRowMax = 72};
  StTpcRSMaker(double eCutOff = 1e-3, const char* name = "TpcRS");
  ~StTpcRSMaker();
  void InitRun(int runnumber = 0);
  void Make(const std::vector<g2t_tpc_hit_st>& g2t_tpc_hit,
                  std::vector<g2t_track_st>& g2t_track,
            const std::vector<g2t_vertex_st>& g2t_vertex, tpcrs::DigiData& digi_data);
  virtual void Finish();
  int         Debug() const {return 1;}
  TF1F* GetPolya(int io = 0)       {return (TF1F*) mPolya[io];}
  TF1F* GetTimeShape0(int io = 0)  {return fgTimeShape0[io];}
  TF1F* GetTimeShape3(int io = 0)  {return fgTimeShape3[io];}
  double GetNoPrimaryClusters(double betaGamma, int charge);
  void Print(Option_t* option = "") const;
  void DigitizeSector(int sector, tpcrs::DigiData& digi_data);
  static int    AsicThresholds(short* ADCs);
  static int    SearchT(const void* elem1, const void** elem2);
  static int    CompareT(const void** elem1, const void** elem2);
  static double shapeEI(double* x, double* par = 0);
  static double shapeEI_I(double* x, double* par = 0);
  static double shapeEI3(double* x, double* par = 0);
  static double shapeEI3_I(double* x, double* par = 0);
  static double fei(double t, double t0, double T);
  static double polya(double* x, double* par);
  SignalSum_t*  GetSignalSum(int sector);
  SignalSum_t*  ResetSignalSum(int sector);
  static double Ec(double* x, double* p); // minimal energy to create an ion pair
  static TF1* fEc(double w = 26.2);           // HEED function to generate Ec
 private:
  static double ShaperFunc(double* x, double* p);
  static double PadResponseFunc(double* x, double* p);
  static double Gatti(double* x, double* p);
  static double InducedCharge(double s, double h, double ra, double Va, double &t0);
  bool TrackSegment2Propagate(g2t_tpc_hit_st* tpc_hitC, const g2t_vertex_st* gver, HitPoint_t &TrackSegmentHits);
  void   GenerateSignal(HitPoint_t &TrackSegmentHits, int sector, int rowMin, int rowMax, double sigmaJitterT, double sigmaJitterX);
  double dEdxCorrection(HitPoint_t &TrackSegmentHits);
  static TF1F*     fgTimeShape3[2];  //!
  static TF1F*     fgTimeShape0[2];   //!
  int    m_Mode;
  SignalSum_t*     m_SignalSum;       //!
  TH1D*    mdNdx;                     //!
  TH1D*    mdNdxL10;                  //!
  TH1D*    mdNdEL10;                  //!
  TF1F*  mShaperResponses[2][24];     //!
  TF1F*  mShaperResponse;             //!
  TF1F*  mChargeFraction[2][24];      //!
  TF1F*  mPadResponseFunction[2][24]; //!
  TF1F*  mPolya[2];                   //!
  TF1F*  mGG;                         //! Gating Grid Transperency
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
  double   msMin, msMax;            //!
  int      mNSplittedHits;          //!
  double xOnWire, yOnWire, zOnWire; //!
  double mGainLocal;                //!
  double QAv;                       //!
  double TotalSignal;               //!
  int pad0;                         //!
  int tbk0;                         //!
  double TotalSignalInCluster;      //!
  double padsdE[kPadMax];           //!
  double tbksdE[kTimeBacketMax];    //!
  double rowsdEH[kRowMax];          //!
  double rowsdE[kRowMax];           //!
  double             mLaserScale;   //!
  const double minSignal;           //!
  const double ElectronRange;       //!
  const double ElectronRangeEnergy; //!
  const double ElectronRangePower;  //!
  const int max_sectors;            //!
  const int max_pads;               //!
  static const int max_timebins = 512;
  double   mCutEle;                 //! cut for delta electrons
};
#endif
