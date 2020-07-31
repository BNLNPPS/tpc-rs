#pragma once

#include "config_structs.h"
#include "struct_containers.h"


struct dE_t {
  float dE;
  float dx;
  float dEdx;
  float dEdxL;
  float dEdxN;
};


struct dEdxCorrection_t {
  dEdxCorrection_t(const char* name = "", const  char* title = "", tpcrs::IConfigStruct* chair = 0, int n = 0) :
    Name(name), Title(title), Chair(chair), nrows(n), dE(0) {}
  const char* Name;
  const char* Title;
  tpcrs::IConfigStruct* Chair;
  int   nrows;
  float dE;
};


class dEdxY2_t;
class StTpcdEdxCorrection
{
 public:
  enum ESector : int {kTpcOuter = 0, kTpcInner = 1, kiTpc = 2};
  enum EOptions : int {
    kUncorrected           =  0,//U
    kAdcCorrection         =  1,//R
    kEdge                  =  2,//E   correction near edge of chamber
    kAdcCorrectionMDF      =  3,//RMDF
    kTpcdCharge            =  4,//D
    kTpcrCharge            =  5,//D
    kTpcCurrentCorrection  =  6,//
    kTpcSecRowB            =  7,//S
    kTpcSecRowC            =  8,//S
    kTpcRowQ               =  9,//
    ktpcPressure           = 10,//P
    ktpcTime               = 11,//t
    kDrift                 = 12,//O
    kMultiplicity          = 13,//M
    kzCorrection           = 14,//Z
    ktpcMethaneIn          = 15,//m
    ktpcGasTemperature     = 16,//T
    ktpcWaterOut           = 17,//W   				7
    kSpaceCharge           = 18,//C   space charge near the wire
    kPhiDirection          = 19,//p   correction wrt local interception angle
    kTanL                  = 20,//p   correction wrt local tan(lambda)
    kdXCorrection          = 21,//X
    kTpcEffectivedX        = 22,//X   Effective pad row height
    kTpcPadTBins           = 23,//d
    kTpcZDC                = 24,//
    kTpcPadMDF             = 25,
    kTpcLast               = 26,//
    kTpcNoAnodeVGainC      = 27,//
    kTpcLengthCorrection   = 28,//
    kTpcLengthCorrectionMDF = 29, //
    kTpcdEdxCor            = 30,//
    kTpcAllCorrections     = 31 //
  };

  struct Corrections {
    enum Bits : unsigned {
      kAll                    = ~(0u),
      kAdcCorrection          = 1u <<  1,
      kEdge                   = 1u <<  2,
      kAdcCorrectionMDF       = 1u <<  3,
      kTpcdCharge             = 1u <<  4,
      kTpcrCharge             = 1u <<  5,
      kTpcCurrentCorrection   = 1u <<  6,
      kTpcSecRowB             = 1u <<  7,
      kTpcSecRowC             = 1u <<  8,
      kTpcRowQ                = 1u <<  9,
      ktpcPressure            = 1u << 10,
      ktpcTime                = 1u << 11,
      kDrift                  = 1u << 12,
      kMultiplicity           = 1u << 13,
      kzCorrection            = 1u << 14,
      ktpcMethaneIn           = 1u << 15,
      ktpcGasTemperature      = 1u << 16,
      ktpcWaterOut            = 1u << 17,
      kSpaceCharge            = 1u << 18,
      kPhiDirection           = 1u << 19,
      kTanL                   = 1u << 20,
      kdXCorrection           = 1u << 21,
      kTpcEffectivedX         = 1u << 22,
      kTpcPadTBins            = 1u << 23,
      kTpcZDC                 = 1u << 24,
      kTpcPadMDF              = 1u << 25,
      kTpcLast                = 1u << 26,
      kTpcNoAnodeVGainC       = 1u << 27,
      kTpcLengthCorrection    = 1u << 28,
      kTpcLengthCorrectionMDF = 1u << 29,
      kTpcdEdxCor             = 1u << 30,
      kTpcAllCorrections      = 1u << 31
    };
  };

  StTpcdEdxCorrection(const tpcrs::Configurator& cfg, int Option = 0, int debug = 0);
  ~StTpcdEdxCorrection();
  int dEdxCorrection(dEdxY2_t &dEdx, bool doIT = true);

  void ReSetCorrections();
  //  St_trigDetSums    *trigDetSums()         {return m_trigDetSums;}

  float           Adc2GeV()              {return mAdc2GeV;}
  void Print(Option_t* opt = "") const;
 private:

  const tpcrs::Configurator& cfg_;
  int                m_Mask;                  //!
  dEdxY2_t*            mdEdx;
  float              mAdc2GeV;                //! Outer/Inner conversion factors from ADC -> GeV
  dEdxCorrection_t     m_Corrections[kTpcAllCorrections];//!
  int                m_Debug;                 //!
};


struct dEdxY2_t
{
  /* U->R->S->P->O->Z->X
     U->R (TpcAdcCorrection) -> P (tpcPressure) ->
     S (TpcSecRowB/TpcSecRowC) ->  O (TpcDriftDistOxygen) ->
     Z (TpcZCorrection) -> X(TpcdXCorrection) */
  int    sector;
  int    row;
  int    channel;
  float  pad;
  int    Npads;
  int    Ntbins;
  float  ZdriftDistance;     // drift distance
  float  ZdriftDistanceO2;   // ZdriftDistance*ppmOxygenIn
  float  ZdriftDistanceO2W;  // ZdriftDistance*ppmOxygenIn*ppmWaterOut
  float  DeltaZ;             // distance to privious cluster
  float  QRatio;             // Ratio to previous cluster Charge
  float  QRatioA;            // Ratio to Sum of all previous cluster Charge
  float  QSumA;              // Sum of all previous cluster Charge
  float  dxC;                // corrected dx which should be used with FitN
  float  xyz[3];  // local
  float  xyzD[3]; // local direction
  float  edge;    // distance to sector edge
  float  PhiR;    // relative phi
  float  resXYZ[3]; // track SectorLocal residual wrt local track
  float  Prob;
  float  zdev;
  float  zP;      // the most probable value from Bichsel
  float  zG;      // global z oh Hit
  float  sigmaP;  // sigma from Bichsel
  float  dCharge; // d_undershoot_Q/Q = ratio of modified - original charge normalized on original charge
  float  rCharge; // d_rounding_Q/Q   = estimated rounding normalized on original charge
  int    lSimulated;
  float  Qcm;     // accumulated charge uC/cm
  float  Crow;    // Current per row;
  float  Zdc;     // ZDC rate from trigger
  float  Weight;  // 1/.sigma^2 of TpcSecRow gas gain correction
  float  adc;     //  adc count from cluster finder
  float  TanL;
  float  Voltage; // Anode Voltage
  float  xpad;    // relative position in pad [-1.0,1.0]
  float  yrow;    // relative position in row [-0.5,0.0] inner, and [0.0,0.5] outer
  float  tpcTime;
  dE_t     C[StTpcdEdxCorrection::kTpcAllCorrections + 1];
  dE_t     F;     //!
  int    npads; // cluster size in pads
  int    ntmbks;// clustre size in time buckets
};
