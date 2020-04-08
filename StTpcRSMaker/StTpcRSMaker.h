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
  Float_t      Sum;
  Short_t      Adc;
  Short_t  TrackId;
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
  virtual              ~StTpcRSMaker();
  virtual void InitRun(int runnumber = 0);
  virtual void Make(const St_g2t_tpc_hit* g2t_tpc_hit, const St_g2t_track* g2t_track, const St_g2t_vertex*  g2t_ver, tpcrs::DigiData& digi_data);
  virtual void Finish();
  Int_t         Debug() const {return 1;}
  TF1F* GetPolya(Int_t io = 0)       {return (TF1F*) mPolya[io];}
  TF1F* GetTimeShape0(Int_t io = 0)  {return fgTimeShape0[io];}
  TF1F* GetTimeShape3(Int_t io = 0)  {return fgTimeShape3[io];}
  Double_t GetNoPrimaryClusters(Double_t betaGamma, Int_t charge);
  virtual void Print(Option_t* option = "") const;
  void DigitizeSector(Int_t sector, tpcrs::DigiData& digi_data);
  static Int_t    AsicThresholds(Short_t* ADCs);
  static Int_t    SearchT(const void* elem1, const void** elem2);
  static Int_t    CompareT(const void** elem1, const void** elem2);
  static Double_t shapeEI(Double_t* x, Double_t* par = 0);
  static Double_t shapeEI_I(Double_t* x, Double_t* par = 0);
  static Double_t shapeEI3(Double_t* x, Double_t* par = 0);
  static Double_t shapeEI3_I(Double_t* x, Double_t* par = 0);
  static Double_t fei(Double_t t, Double_t t0, Double_t T);
  static Double_t polya(Double_t* x, Double_t* par);
  SignalSum_t*  GetSignalSum(Int_t sector);
  SignalSum_t*  ResetSignalSum(Int_t sector);
  static Double_t Ec(Double_t* x, Double_t* p); // minimal energy to create an ion pair
  static TF1* fEc(Double_t w = 26.2);           // HEED function to generate Ec
 private:
  static Double_t ShaperFunc(Double_t* x, Double_t* p);
  static Double_t PadResponseFunc(Double_t* x, Double_t* p);
  static Double_t Gatti(Double_t* x, Double_t* p);
  static Double_t InducedCharge(Double_t s, Double_t h, Double_t ra, Double_t Va, Double_t &t0);
  Bool_t TrackSegment2Propagate(g2t_tpc_hit_st* tpc_hitC, g2t_vertex_st* gver, HitPoint_t &TrackSegmentHits);
  void   GenerateSignal(HitPoint_t &TrackSegmentHits, Int_t sector, Int_t rowMin, Int_t rowMax, Double_t sigmaJitterT, Double_t sigmaJitterX);
  Double_t dEdxCorrection(HitPoint_t &TrackSegmentHits);
  static TF1F*     fgTimeShape3[2];  //!
  static TF1F*     fgTimeShape0[2];   //!
  Int_t    m_Mode;
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
  Double_t InnerAlphaVariation[24];   //!
  Double_t OuterAlphaVariation[24];   //!
  Altro* mAltro;                      //!
  Int_t    numberOfInnerSectorAnodeWires; //!
  Double_t firstInnerSectorAnodeWire; //!
  Double_t lastInnerSectorAnodeWire;  //!
  Int_t    numberOfOuterSectorAnodeWires; //!
  Double_t firstOuterSectorAnodeWire; //!
  Double_t lastOuterSectorAnodeWire;  //!
  Double_t anodeWirePitch;            //!
  Double_t anodeWireRadius;           //!
  Double_t innerSectorAnodeVoltage[24];//!
  Double_t outerSectorAnodeVoltage[24];//!
  Double_t    mLocalYDirectionCoupling[2][24][7]; //!
  Double_t   msMin, msMax;            //!
  Int_t      mNSplittedHits;          //!
  Double_t xOnWire, yOnWire, zOnWire; //!
  Double_t mGainLocal;                //!
  Double_t QAv;                       //!
  Double_t TotalSignal;               //!
  Int_t pad0;                         //!
  Int_t tbk0;                         //!
  Double_t TotalSignalInCluster;      //!
  Double_t padsdE[kPadMax];           //!
  Double_t tbksdE[kTimeBacketMax];    //!
  Double_t rowsdEH[kRowMax];          //!
  Double_t rowsdE[kRowMax];           //!
  Double_t             mLaserScale;   //!
  const Double_t minSignal;           //!
  const Double_t ElectronRange;       //!
  const Double_t ElectronRangeEnergy; //!
  const Double_t ElectronRangePower;  //!
  const Int_t max_sectors;            //!
  const Int_t max_pads;               //!
  static const Int_t max_timebins = 512;
  Double_t   mCutEle;                 //! cut for delta electrons
};
#endif
