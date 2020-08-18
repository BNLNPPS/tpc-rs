#include <cassert>
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

#include "logger.h"
#include "tpcrs/detail/struct_containers.h"

#define MakeChairInstance(STRUCT,PATH) \
template<> std::string tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_ ## STRUCT ## C , STRUCT>::name(# PATH);

#define MakeChairInstance2(STRUCT,CLASS,PATH) \
template<> std::string tpcrs::ConfigStruct<St_ ## STRUCT ## C, CLASS , STRUCT>::name(# PATH);

MakeChairInstance(TpcEffectivedX,Calibrations/tpc/TpcEffectivedX);
MakeChairInstance(tpcGridLeak,Calibrations/tpc/tpcGridLeak);
MakeChairInstance2(TpcSecRowCor,St_TpcSecRowBC,Calibrations/tpc/TpcSecRowB);
MakeChairInstance2(TpcSecRowCor,St_TpcSecRowCC,Calibrations/tpc/TpcSecRowC);
MakeChairInstance(tpcSCGL,Calibrations/tpc/tpcSCGL);


double St_tpcCorrectionC::CalcCorrection(int i, double x, double z, int NparMax) const {
  tpcCorrection *cor =  Struct(i);
  return SumSeries(cor, x, z, NparMax);
}


double St_tpcCorrectionC::SumSeries(tpcCorrection *cor,  double x, double z, int NparMax) const {
  double Sum = 0;
  if (! cor) return Sum;
  int N = std::abs(cor->npar)%100;
  if (N == 0) return Sum;
  if (NparMax > 0) N = NparMax;
  static double T0, T1, T2;
  // parameterization variable
  double X = x;
  if (cor->npar  < 0) X = std::exp(x);
  else {
    switch  (cor->type) {
    case 10:// ADC correction offset + poly for ADC
    case 11:// ADC correction offset + poly for log(ADC) and |Z|
    case 12:// ADC correction offset + poly for log(ADC) and TanL
      X = std::log(x);      break;
    case 1: // Tchebyshev [-1,1]
      if (cor->min < cor->max)   X = -1 + 2*std::max(0.,std::min(1.,(X - cor->min)/( cor->max - cor->min)));
      break;
    case 2: // Shifted TChebyshev [0,1]
      if (cor->min < cor->max)   X = std::max(0.,std::min(1.,(X - cor->min)/( cor->max - cor->min)));
      break;
    case 3:
      if (std::abs(x) >= 1) X = 0;
      else                    X = std::log(1. - std::abs(x));
      break;
    case 4:
      if (std::abs(x) >= 1) X = 0;
      else                    X = std::copysign(std::log(1. - std::abs(x)),x);
      break;
    case 5:
      if (x < 1e-7) X = -16.118;
      else          X = std::log(x);
      break;
    default:      X = x;    break;
    }
  }
  if (cor->type != 1 && cor->type != 2 &&
      cor->min < cor->max) {
    if (X < cor->min) X = cor->min;
    if (X > cor->max) X = cor->max;
  }
  static TF1 *f1000 = 0, *f1100 = 0, *f1200 = 0, *f1300 = 0;
  static TF1 *f2000 = 0, *f2100 = 0, *f2200 = 0, *f2300 = 0;
  TF1 *f = 0;
  switch (cor->type) {
  case 1: // Tchebyshev [-1,1]
    T0 = 1;
    Sum = cor->a[0]*T0;
    if (N == 1) break;
    T1 = X;
    Sum += cor->a[1]*T1;
    for (int n = 2; n <= N; n++) {
      T2 = 2*X*T1 - T0;
      Sum += cor->a[n]*T2;
      T0 = T1;
      T1 = T2;
    }
    break;
  case 2: // Shifted TChebyshev [0,1]
    T0 = 1;
    Sum = cor->a[0]*T0;
    if (N == 1) break;
    T1 = 2*X - 1;
    Sum += cor->a[1]*T1;
    for (int n = 2; n <= N; n++) {
      T2 = 2*(2*X - 1)*T1 - T0;
      Sum += cor->a[n]*T2;
      T0 = T1;
      T1 = T2;
    }
    break;
  case 10: // ADC correction offset + poly for ADC
    Sum = cor->a[N-1];
    for (int n = N-2; n>=0; n--) Sum = X*Sum + cor->a[n];
    Sum += std::log(1. + cor->OffSet/x);
    Sum  = std::exp(Sum);
    Sum *= x;
    break;
  case 11: // ADC correction offset + poly for log(ADC) and |Z|
    Sum = cor->a[1] + z*cor->a[2] + z*X*cor->a[3] + std::exp(X*(cor->a[4] + X*cor->a[5]) + cor->a[6]);
    Sum *= std::exp(-cor->a[0]);
    break;
  case 12: // ADC correction offset + poly for log(ADC) and TanL
    Sum = cor->a[1] + z*cor->a[2] + z*z*cor->a[3] + std::exp(X*(cor->a[4] + X*cor->a[5]) + cor->a[6]);
    Sum *= std::exp(-cor->a[0]);
    break;
  case 1000:
  case 1100:
  case 1200:
  case 1300:
    if (cor->type == 1000) {
      if (! f1000) f1000 = new TF1("f1000","gaus+pol0(3)");
      f = f1000;
    } else if (cor->type == 1100) {
      if (! f1100) f1100 = new TF1("f1100","gaus+pol1(3)");
      f = f1100;
    } else if (cor->type == 1200) {
      if (! f1200) f1200 = new TF1("f1200","gaus+pol2(3)");
      f = f1200;
    } else if (cor->type == 1300) {
      if (! f1300) f1300 = new TF1("f1300","gaus+pol3(3)");
      f = f1300;
    }
    assert(f);
    f->SetParameters(cor->a);
    Sum = f->Eval(X);
    break;
  case 2000:
  case 2100:
  case 2200:
  case 2300:
    if (cor->type == 2000) {
      if (! f2000) f2000 = new TF1("f2000","expo+pol0(2)");
      f = f2000;
    } else if (cor->type == 2100) {
      if (! f2100) f2100 = new TF1("f2100","expo+pol1(2)");
      f = f2100;
    } else if (cor->type == 2200) {
      if (! f2200) f2200 = new TF1("f2200","expo+pol2(2)");
      f = f2200;
    } else if (cor->type == 2300) {
      if (! f2300) f2300 = new TF1("f2300","expo+pol3(2)");
      f = f2300;
    }
    assert(f);
    f->SetParameters(cor->a);
    Sum = f->Eval(X);
    break;
  default: // polynomials
    Sum = cor->a[N-1];
    for (int n = N-2; n>=0; n--) Sum = X*Sum + cor->a[n];
    break;
  }
  return Sum;
}

MakeChairInstance2(tpcCorrection,St_TpcRowQC,Calibrations/tpc/TpcRowQ);
MakeChairInstance2(tpcCorrection,St_TpcDriftDistOxygenC,Calibrations/tpc/TpcDriftDistOxygen);
MakeChairInstance2(tpcCorrection,St_TpcMultiplicityC,Calibrations/tpc/TpcMultiplicity);
MakeChairInstance2(tpcCorrection,St_TpcZCorrectionBC,Calibrations/tpc/TpcZCorrectionB);
MakeChairInstance2(tpcCorrection,St_TpcdXCorrectionBC,Calibrations/tpc/TpcdXCorrectionB);
MakeChairInstance2(tpcCorrection,St_tpcPressureBC,Calibrations/tpc/tpcPressureB);
MakeChairInstance2(tpcCorrection,St_TpcEdgeC,Calibrations/tpc/TpcEdge);
MakeChairInstance2(tpcCorrection,St_TpcAdcCorrectionBC,Calibrations/tpc/TpcAdcCorrectionB);
MakeChairInstance2(MDFCorrection,St_TpcAdcCorrectionMDF,Calibrations/tpc/TpcAdcCorrectionMDF);
MakeChairInstance2(tpcCorrection,St_tpcMethaneInC,Calibrations/tpc/tpcMethaneIn);
MakeChairInstance2(tpcCorrection,St_tpcGasTemperatureC,Calibrations/tpc/tpcGasTemperature);
MakeChairInstance2(tpcCorrection,St_tpcWaterOutC,Calibrations/tpc/tpcWaterOut);
MakeChairInstance2(tpcCorrection,St_tpcTimeDependenceC,Calibrations/tpc/tpcTimeDependence);
MakeChairInstance2(tpcCorrection,St_TpcdChargeC,Calibrations/tpc/TpcdCharge);
MakeChairInstance2(tpcCorrection,St_TpcrChargeC,Calibrations/tpc/TpcrCharge);
MakeChairInstance2(tpcCorrection,St_TpcTanLC,Calibrations/tpc/TpcTanL);
MakeChairInstance2(tpcCorrection,St_TpcCurrentCorrectionC,Calibrations/tpc/TpcCurrentCorrectionX);
MakeChairInstance2(tpcCorrection,St_TpcZDCC,Calibrations/tpc/TpcZDC);
MakeChairInstance2(tpcCorrection,St_TpcSpaceChargeC,Calibrations/tpc/TpcSpaceCharge);
MakeChairInstance2(tpcCorrection,St_TpcPhiDirectionC,Calibrations/tpc/TpcPhiDirection);
MakeChairInstance2(tpcCorrection,St_TpcLengthCorrectionBC,Calibrations/tpc/TpcLengthCorrectionB);
MakeChairInstance2(MDFCorrection,St_TpcLengthCorrectionMDF,Calibrations/tpc/TpcLengthCorrectionMDF);
MakeChairInstance2(MDFCorrection,St_TpcPadCorrectionMDF,Calibrations/tpc/TpcPadCorrectionMDF);
St_MDFCorrectionC *St_MDFCorrectionC::fgMDFCorrectionC = 0;


St_MDFCorrectionC::St_MDFCorrectionC() : fFunc(0) {
}


St_MDFCorrectionC::~St_MDFCorrectionC() {
  unsigned int N = GetNRows();
  for (unsigned int i = 0; i < N; i++) {SafeDelete(fFunc[i]);}
  delete [] fFunc;
}


double St_MDFCorrectionC::MDFunc(double *x, double *p) {
  // Evaluate parameterization at point x. Optional argument coeff is
  // a vector of coefficients for the parameterisation, NCoefficients
  // elements long.
  assert(x);
  unsigned int k = p[0];
  assert(k >= 0 && k < fgMDFCorrectionC->GetNRows());
  double returnValue = fgMDFCorrectionC->DMean(k);
  double term        = 0;
  unsigned char    i, j;
  for (i = 0; i < fgMDFCorrectionC->NCoefficients(k); i++) {
    // Evaluate the ith term in the expansion
    term = fgMDFCorrectionC->Coefficients(k)[i];
    for (j = 0; j < fgMDFCorrectionC->NVariables(k); j++) {
      // Evaluate the factor (polynomial) in the j-th variable.
      int    p  =  fgMDFCorrectionC->Powers(k)[i * fgMDFCorrectionC->NVariables(k) + j];
      double y  =  1 + 2. / (fgMDFCorrectionC->XMax(k)[j] - fgMDFCorrectionC->XMin(k)[j])
	* (x[j] - fgMDFCorrectionC->XMax(k)[j]);
      term        *= fgMDFCorrectionC->EvalFactor(k,p,y);
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}



double St_MDFCorrectionC::Eval(int k, double x0, double x1) const {
  double x[2] = {x0, x1};
  return Eval(k,x);
}


double St_MDFCorrectionC::Eval(int k, double *x) const {
  // Evaluate parameterization at point x. Optional argument coeff is
  // a vector of coefficients for the parameterisation, NCoefficients
  // elements long.
  assert(x);
  if (! fFunc[k]) {
    fgMDFCorrectionC = (St_MDFCorrectionC *) this;
    if (NVariables(k) <= 0) {
      return 0;
    } else if (NVariables(k) == 1) {
      fFunc[k] = new TF1(Form("%s_%i",GetName().c_str(),k),St_MDFCorrectionC::MDFunc,
			 XMin(k)[0],XMax(k)[0],1);
      fFunc[k]->SetParameter(0,k);
      fFunc[k]->Save(XMin(k)[0],XMax(k)[0],0,0,0,0);
    } else if (NVariables(k) == 2) {
      fFunc[k] = new TF2(Form("%s_%i",GetName().c_str(),k),St_MDFCorrectionC::MDFunc,
			 XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],1);
      fFunc[k]->SetParameter(0,k);
      ((TF2 *) fFunc[k])->Save(XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],0,0);
    } else if (NVariables(k) == 3) {
      fFunc[k] = new TF3(Form("%s_%i",GetName().c_str(),k),St_MDFCorrectionC::MDFunc,
			 XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],XMin(k)[2],XMax(k)[2],1);
      fFunc[k]->SetParameter(0,k);
      ((TF3 *) fFunc[k])->Save(XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],XMin(k)[2],XMax(k)[2]);
    }
  }
  double xx[3];
  for (int v = 0; v < NVariables(k); v++) {
    xx[v] = std::max(XMin(k)[v], std::min(XMin(k)[v]+0.999*(XMax(k)[v]-XMin(k)[v]), x[v]));
  }
  return fFunc[k]->GetSave(xx);
}


double St_MDFCorrectionC::EvalError(int k, double *x) const {
  // Evaluate parameterization error at point x. Optional argument coeff is
  // a vector of coefficients for the parameterisation, NCoefficients(k)
  // elements long.
  assert(x);
  double returnValue = 0;
  double term        = 0;
  unsigned char    i, j;
  for (i = 0; i < NCoefficients(k); i++) {
    // Evaluate the ith term in the expansion
    term = CoefficientsRMS(k)[i];
    for (j = 0; j < NVariables(k); j++) {
      // Evaluate the factor (polynomial) in the j-th variable.
      int    p  =  Powers(k)[i * NVariables(k) + j];
      double y  =  1 + 2. / (XMax(k)[j] - XMin(k)[j]) * (x[j] - XMax(k)[j]);
      term     *= EvalFactor(p,y);
    }
    // Add this term to the final result
    returnValue += term*term;
  }
  return std::sqrt(returnValue);
}


double St_MDFCorrectionC::EvalFactor(int k, int p, double x) const {
  // Evaluate function with power p at variable value x
  int    i   = 0;
  double p1  = 1;
  double p2  = 0;
  double p3  = 0;
  double r   = 0;

  switch(p) {
  case 1:
    r = 1;
    break;
  case 2:
    r =  x;
    break;
  default:
    p2 = x;
    for (i = 3; i <= p; i++) {
      p3 = p2 * x;
      if (PolyType(k) == kLegendre)
	p3 = ((2 * i - 3) * p2 * x - (i - 2) * p1) / (i - 1);
      else if (PolyType(k) == kChebyshev)
	p3 = 2 * x * p2 - p1;
      p1 = p2;
      p2 = p3;
    }
    r = p3;
  }
  return r;
}

MakeChairInstance(tpcAnodeHV,Calibrations/tpc/tpcAnodeHV);

void  St_tpcAnodeHVC::sockets(int sector, int padrow, int &e1, int &e2, float &f2) {
  e1 = (sector-1)*19;
  e2 = e1;
  f2 = 0;
  // sector=1..24 , padrow=1..45
  // f2 represents signal couplings from neighboring HV sections
  // see: http://www.star.bnl.gov/public/tpc/hard/signals/signal_division.html
    switch (padrow) {
    case  1: e1+= 1; e2+= 2; f2 = 0.00197; break;
    case  2: e1+= 2; break;
    case  3: e1+= 3; e2+= 2; f2 = 0.04547; break;
    case  4: e1+= 3; break;
    case  5: e1+= 4; break;
    case  6: e1+= 4; e2+= 5; f2 = 0.00007; break;
    case  7: e1+= 5; break;
    case  8: e1+= 6; e2+= 5; f2 = 0.04547; break;
    case  9: e1+= 6; break;
    case 10: e1+= 7; break;
    case 11: e1+= 8; e2+= 7; f2 = 0.33523; break;
    case 12: e1+= 8; break;
    case 13: e1+=17; break;
    case 14: e1+= 9; e2+=18; f2 = 0.00312; break;
    case 15:
    case 16: e1+= 9; break;
    case 17: e1+= 9; e2+=10; f2 = 0.40250; break;
    case 18:
    case 19:
    case 20: e1+=10; break;
    case 21: e1+=10; e2+=11; f2 = 0.40250; break;
    case 22:
    case 23:
    case 24: e1+=11; break;
    case 25: e1+=11; e2+=12; f2 = 0.40250; break;
    case 26:
    case 27:
    case 28: e1+=12; break;
    case 29: e1+=12; e2+=13; f2 = 0.40250; break;
    case 30:
    case 31:
    case 32: e1+=13; break;
    case 33: e1+=13; e2+=14; f2 = 0.40250; break;
    case 34:
    case 35:
    case 36: e1+=14; break;
    case 37: e1+=14; e2+=15; f2 = 0.40250; break;
    case 38:
    case 39:
    case 40: e1+=15; break;
    case 41: e1+=15; e2+=16; f2 = 0.40250; break;
    case 42:
    case 43:
    case 44: e1+=16; break;
    case 45: e1+=16; e2+=19; f2 = 0.40250; break;
    default: e1 = 0; e2 = 0; f2 = 0;
    }
}


MakeChairInstance(TpcAvgPowerSupply,Calibrations/tpc/TpcAvgPowerSupply);


float St_TpcAvgPowerSupplyC::voltagePadrow(int sector, int padrow) const {
  int e1 = 0, e2 = 0;
  float f2 = 0;
  St_tpcAnodeHVC::sockets(sector, padrow, e1, e2, f2);
  if (e1==0) return -99;
  int ch1 = St_TpcAvgCurrentC::ChannelFromSocket((e1-1)%19+1);
  float v1=Voltage()[8*(sector-1)+ch1-1] ;
  if (f2==0) return v1;
  int ch2 = St_TpcAvgCurrentC::ChannelFromSocket((e2-1)%19 + 1);
  if (ch1 == ch2) return v1;
  float v2=Voltage()[8*(sector-1)+ch2-1] ;
  if (v2==v1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  float B = tpcrs::IsInner(padrow, cfg_) ? 13.05e-3 : 10.26e-3;
  float v_eff = std::log((1.0-f2)*std::exp(B*v1) + f2*std::exp(B*v2)) / B;
  return v_eff;
}


float St_TpcAvgPowerSupplyC::AcChargeL(int sector, int channel) {
  // L Inner = 190222, Outer = 347303
  static float Length[8] = {
    1307.59, //   Channel 1
    1650.57, //   Channel 2
    1993.54, //   Channel 3
    2974.24, //   Channel 4
    3324.59, //   Channel 5
    3202.42, //   Channel 6
    3545.4 , //   Channel 7
    4398.53};//   Channel 8

  return AcCharge(sector,channel)/Length[channel-1];
}

MakeChairInstance2(tpcCorrection,St_tpcGainCorrectionC,Calibrations/tpc/tpcGainCorrection);
MakeChairInstance(TpcAvgCurrent,Calibrations/tpc/TpcAvgCurrent);


int St_TpcAvgCurrentC::ChannelFromSocket(int socket) {
  int channel = -1;
  switch (socket) {
  case 1:
  case 2 : channel = 1; break;
  case 3:
  case 4:  channel = 2; break;
  case 5:
  case 6:  channel = 3; break;
  case 7:
  case 8:
  case 17: channel = 4; break;
  case 9:
  case 10:
  case 18: channel = 5; break;
  case 11:
  case 12: channel = 6; break;
  case 13:
  case 14: channel = 7; break;
  case 15:
  case 16:
  case 19: channel = 8; break;
  default:              break;
  }
  return channel;
}


float St_TpcAvgCurrentC::AcChargeL(int sector, int channel) {
  if (! cfg_.C<St_TpcAvgPowerSupplyC>().is_missing) {
    return cfg_.C<St_TpcAvgPowerSupplyC>().AcChargeL(sector,channel);
  }
  static float Length[8] = {
    1307.59, //   Channel 1
    1650.57, //   Channel 2
    1993.54, //   Channel 3
    2974.24, //   Channel 4
    3324.59, //   Channel 5
    3202.42, //   Channel 6
    3545.4 , //   Channel 7
    4398.53};//   Channel 8
  return AcCharge(sector,channel)/Length[channel-1];
}


float St_TpcAvgCurrentC::AvCurrent(int sector, int channel) {
  if (! cfg_.C<St_TpcAvgPowerSupplyC>().is_missing) {
    return cfg_.C<St_TpcAvgPowerSupplyC>().AvCurrent(sector,channel);
  }
  return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
    Struct()->AvCurrent[8*(sector-1)+channel-1] :     0;
}


float St_TpcAvgCurrentC::AcCharge(int sector, int channel) {
  if (! cfg_.C<St_TpcAvgPowerSupplyC>().is_missing) {
    return cfg_.C<St_TpcAvgPowerSupplyC>().AcCharge(sector,channel);
  }
  return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
    Struct()->AcCharge[8*(sector-1)+channel-1] :     0;
}
MakeChairInstance(tpcRDOMap,Calibrations/tpc/tpcRDOMap);


int St_tpcRDOMapC::rdo(int padrow, int pad) const {
  int rdo = 0;
  int N = nrows(0);
  for (int i = 0; i < N; i++) {
    if (padrow != row(i)) continue;
    if (pad < padMin(i) || pad > padMax(i)) continue;
    rdo = rdoI(i);
    
    break;
  }
  return rdo;
}

MakeChairInstance2(spaceChargeCor,St_spaceChargeCorR1C,Calibrations/rich/spaceChargeCor);
MakeChairInstance2(spaceChargeCor,St_spaceChargeCorR2C,Calibrations/rich/spaceChargeCorR2);
MakeChairInstance(tpcFieldCage,Geometry/tpc/tpcFieldCage);
MakeChairInstance(tpcPadPlanes,Geometry/tpc/tpcPadPlanes);
MakeChairInstance(tpcGlobalPosition,Geometry/tpc/tpcGlobalPosition);
MakeChairInstance(tpcFieldCageShort,Geometry/tpc/tpcFieldCageShort);
MakeChairInstance(tpcHVPlanes,Geometry/tpc/tpcHVPlanes);
MakeChairInstance2(Survey,StTpcInnerSectorPosition,Geometry/tpc/TpcInnerSectorPosition);
MakeChairInstance2(Survey,StTpcOuterSectorPosition,Geometry/tpc/TpcOuterSectorPosition);
MakeChairInstance2(Survey,StTpcSuperSectorPosition,Geometry/tpc/TpcSuperSectorPosition);
MakeChairInstance2(Survey,StTpcHalfPosition,Geometry/tpc/TpcHalfPosition);
MakeChairInstance2(Survey,StTpcPosition,Geometry/tpc/TpcPosition);

St_SurveyC::St_SurveyC() : fRotations()  { }

void St_SurveyC::Initialize()
  {
    unsigned int N = GetNRows();
    for (unsigned int i = 0; i < N; i++) {
      fRotations.push_back(new TGeoHMatrix);
      TGeoHMatrix &rot = *fRotations[i];
      if (N == 1) rot.SetName(GetName().c_str());
      else        rot.SetName(Form("%s_%i",GetName().c_str(),i+1));
      rot.SetRotation(Rotation(i));
      rot.SetTranslation(Translation(i));
      assert(TMath::Abs(rot.Determinant())-1 < 1.e-3);
    }
  }


St_SurveyC::~St_SurveyC() {
    for (unsigned int i = 0; i < GetNRows(); i++) {
      SafeDelete(fRotations[i]);
    }
}


const TGeoHMatrix &St_SurveyC::GetMatrix(int i) {
  assert(fRotations[i]);
  assert(std::abs(fRotations[i]->Determinant())-1 < 1.e-3);
  return *fRotations[i];
}


double St_spaceChargeCorC::getSpaceChargeCoulombs(const tpcrs::Configurator& cfg)
  {
    const trigDetSums& scalers = cfg.S<trigDetSums>();

    double coulombs = 0;
    bool use_powers = true;

    for (int row=0;row< (int) GetNRows();row++) {
      double mult = 0;
      switch ((int) getSpaceChargeDetector(row)) {
        case (0) : mult = scalers.mult; break; // vpdx as of 2007-12-19
        case (1) : mult = scalers.bbcX; break;
        case (2) : mult = scalers.zdcX; break;
        case (3) : mult = scalers.zdcEast + scalers.zdcWest; break;
        case (4) : mult = scalers.bbcEast + scalers.bbcWest; break;
        case (5) : mult = scalers.zdcEast; break;
        case (6) : mult = scalers.zdcWest; break;
        case (7) : mult = scalers.bbcEast; break;
        case (8) : mult = scalers.bbcWest; break;
        case (9) : mult = scalers.bbcYellowBkg; break;
        case (10): mult = scalers.bbcBlueBkg; break;
        case (11): mult = scalers.pvpdEast; break;
        case (12): mult = scalers.pvpdWest; break;
        case (13): mult = scalers.ctbTOFp; break; // zdcx-no-killer as of 2011
        case (14): mult = scalers.ctbEast; break; // zdce-no-killer as of 2011
        case (15): mult = scalers.ctbWest; break; // zdcw-no-killer as of 2011

        default  : mult = 0.;
      }
      if (mult < 0) {
        is_missing = true;
        return 0; // Unphysical scaler rates will be uncorrected
      } else
        is_missing = false;

      double saturation = getSpaceChargeSatRate(row);
      double correction = getSpaceChargeCorrection(cfg.S<MagFactor>().ScaleFactor,row);
      double factor     = getSpaceChargeFactor(row);
      double offset     = getSpaceChargeOffset(row);
      double intens = (mult < saturation) ? mult : saturation;

      if (use_powers) coulombs += ::pow(intens-offset,factor) * correction ;
      else coulombs += factor * (intens-offset) * correction ;
    }
    return coulombs;
  }


namespace tpcrs {

int ChannelFromRow(int row)
{
    if (row <  3) return 1;
    if (row <  7) return 2;
    if (row < 10) return 3;
    if (row < 14) return 4;
    if (row < 22) return 5;
    if (row < 30) return 6;
    if (row < 38) return 7;
    return 8;
}

float VoltagePadrow(int sector, int row, const Configurator& cfg)
{
  if (! cfg.C<St_TpcAvgPowerSupplyC>().is_missing) {
    return cfg.C<St_TpcAvgPowerSupplyC>().voltagePadrow(sector,row);
  }
  int e1 = 0, e2 = 0;
  float f2 = 0;
  St_tpcAnodeHVC::sockets(sector, row, e1, e2, f2);
  if (e1==0) return -99;
  float v1=cfg.S<tpcAnodeHVavg>(e1-1).voltage;
  if (f2==0) return v1;
  float v2=cfg.S<tpcAnodeHVavg>(e2-1).voltage;
  if (v2==v1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  float B = tpcrs::IsInner(row, cfg) ? 13.05e-3 : 10.26e-3;
  float v_eff = std::log((1.0-f2)*std::exp(B*v1) + f2*std::exp(B*v2)) / B;
  return v_eff;
}

float GainCorrection(int sector, int row, const tpcrs::Configurator& cfg)
{
  int l = 0;
  double V_nominal = 1390;
  float V = 0;
  float gain = 0;
  if (IsInner(row, cfg)) {l = 1; V_nominal = 1170;}
  const St_tpcGainCorrectionC& gC = cfg.C<St_tpcGainCorrectionC>();

  int NRows = gC.GetNRows();
  if (l >= NRows) return gain;

  V = VoltagePadrow(sector, row, cfg);
  if (V > 0) {
    double v = V - V_nominal;
    if (v < gC.min(l) || v > gC.max(l))
      return 0;

    if (gC.min(l) < -450) {
      // if range was expanded below 150 V then use only the linear approximation
      gain  = std::exp(gC.CalcCorrection(l,v, 0., 2));
    } else {
      gain  = std::exp(gC.CalcCorrection(l,v));
    }
  }
  return gain;
}


float DriftVelocity(int sector, const Configurator& cfg)
{
  const tpcDriftVelocity& dv = cfg.S<tpcDriftVelocity>();
  float drift_velocity;

  if (DetectorSide(sector, cfg) == TPC::Half::first)
    drift_velocity = dv.laserDriftVelocityWest > 0 ? dv.laserDriftVelocityWest : dv.cathodeDriftVelocityWest;
  else
    drift_velocity = dv.laserDriftVelocityEast > 0 ? dv.laserDriftVelocityEast : dv.cathodeDriftVelocityEast;

  return 1e6 * drift_velocity;
}


TPC::Half DetectorSide(int sector, const Configurator& cfg)
{
  return sector <= 12 ? TPC::Half::first : TPC::Half::second;
}


bool IsInner(int row, const Configurator& cfg)
{
  return row <= cfg.S<tpcPadPlanes>().innerPadRows;
}


double RadialDistanceAtRow(int row, const Configurator& cfg)
{
  if (IsInner(row, cfg)) {
    return cfg.S<tpcPadPlanes>().innerRowRadii[row - 1];
  } else {
    int n_inner_rows = cfg.S<tpcPadPlanes>().innerPadRows;
    return cfg.S<tpcPadPlanes>().outerRowRadii[row - 1 - n_inner_rows];
  }
}

}
