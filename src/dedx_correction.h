#pragma once

#include "config_structs.h"
#include "struct_containers.h"

class dEdxY2_t;

class StTpcdEdxCorrection
{
 public:

  enum EOptions : int {
    kUncorrected            =  0, //U
    kAdcCorrection          =  1, //R
    kEdge                   =  2, //E   correction near edge of chamber
    kAdcCorrectionMDF       =  3, //RMDF
    kTpcdCharge             =  4, //D
    kTpcrCharge             =  5, //D
    kTpcCurrentCorrection   =  6,
    kTpcSecRowB             =  7, //S
    kTpcSecRowC             =  8, //S
    kTpcRowQ                =  9,
    ktpcPressure            = 10, //P
    ktpcTime                = 11, //t
    kDrift                  = 12, //O
    kMultiplicity           = 13, //M
    kzCorrection            = 14, //Z
    ktpcMethaneIn           = 15, //m
    ktpcGasTemperature      = 16, //T
    ktpcWaterOut            = 17, //W   7
    kSpaceCharge            = 18, //C   space charge near the wire
    kPhiDirection           = 19, //p   correction wrt local interception angle
    kTanL                   = 20, //p   correction wrt local tan(lambda)
    kdXCorrection           = 21, //X
    kTpcEffectivedX         = 22, //X   Effective pad row height
    kTpcPadTBins            = 23, //d
    kTpcZDC                 = 24,
    kTpcPadMDF              = 25,
    kTpcLast                = 26,
    kTpcNoAnodeVGainC       = 27,
    kTpcLengthCorrection    = 28,
    kTpcLengthCorrectionMDF = 29,
    kTpcdEdxCor             = 30,
    kTpcAllCorrections      = 31 
  };

  struct Corrections {
    enum Bits : unsigned {
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
      kTpcZDC                 = 1u << 24,
      kTpcPadMDF              = 1u << 25,
      kTpcNoAnodeVGainC       = 1u << 27,
      kTpcLengthCorrection    = 1u << 28,
      kTpcLengthCorrectionMDF = 1u << 29
    };
  };

  StTpcdEdxCorrection(const tpcrs::Configurator& cfg,
    unsigned options =
      Corrections::kEdge                   |
      Corrections::kTpcdCharge             |
      Corrections::kTpcrCharge             |
      Corrections::kTpcCurrentCorrection   |
      Corrections::kTpcSecRowB             |
      Corrections::kTpcSecRowC             |
      Corrections::kTpcRowQ                |
      Corrections::ktpcPressure            |
      Corrections::ktpcTime                |
      Corrections::kDrift                  |
      Corrections::kMultiplicity           |
      Corrections::kzCorrection            |
      Corrections::ktpcMethaneIn           |
      Corrections::ktpcGasTemperature      |
      Corrections::ktpcWaterOut            |
      Corrections::kSpaceCharge            |
      Corrections::kPhiDirection           |
      Corrections::kTanL                   |
      Corrections::kTpcEffectivedX         |
      Corrections::kTpcZDC                 |
      Corrections::kTpcPadMDF              |
      Corrections::kTpcLengthCorrection    |
      Corrections::kTpcLengthCorrectionMDF
  );

  int dEdxCorrection(int sector, int row, dEdxY2_t &dEdx);

 private:

  enum ESector : int {kTpcOuter = 0, kTpcInner = 1};

  struct dEdxCorrection_t
  {
    std::string name;
    std::string title;
    tpcrs::IConfigStruct* Chair;
    int   nrows;
    float dE;
  };

  const tpcrs::Configurator& cfg_;
  const int         options_;
  dEdxCorrection_t  corrections_[kTpcAllCorrections];
};


struct dEdxY2_t
{
  struct dE_t {
    float dE;
    float dx;
  };

  /* U->R->S->P->O->Z->X
     U->R (TpcAdcCorrection) -> P (tpcPressure) ->
     S (TpcSecRowB/TpcSecRowC) ->  O (TpcDriftDistOxygen) ->
     Z (TpcZCorrection) -> X(TpcdXCorrection) */
  float  ZdriftDistance;     // drift distance
  float  DeltaZ;             // distance to privious cluster
  float  QRatio;             // Ratio to previous cluster Charge
  float  QRatioA;            // Ratio to Sum of all previous cluster Charge
  float  dxC;                // corrected dx which should be used with FitN
  float  xyz[3];  // local
  float  xyzD[3]; // local direction
  float  edge;    // distance to sector edge
  float  PhiR;    // relative phi
  float  zG;      // global z oh Hit
  float  dCharge; // d_undershoot_Q/Q = ratio of modified - original charge normalized on original charge
  float  rCharge; // d_rounding_Q/Q   = estimated rounding normalized on original charge
  float  Zdc;     // ZDC rate from trigger
  float  adc;     //  adc count from cluster finder
  float  TanL;
  float  xpad;    // Used in kTpcPadMDF: relative position in pad [-1.0,1.0]
  float  yrow;    // Used in kTpcPadMDF: relative position in row [-0.5,0.0] inner, and [0.0,0.5] outer
  float  tpcTime;
  dE_t   C[StTpcdEdxCorrection::kTpcAllCorrections + 1];
  dE_t   F;     // Final overall correction
  int    npads; // Used in kAdcCorrectionMDF cluster size in pads
  int    ntmbks;// Used in kAdcCorrectionMDF clustre size in time buckets
};
