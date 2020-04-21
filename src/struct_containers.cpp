#include <cassert>
#include "TEnv.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "tpcrs/logger.h"

#define MakeChairInstance(STRUCT,PATH) \
template<> std::string tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_ ## STRUCT ## C , STRUCT ## _st>::name(# PATH);

#define MakeChairInstance2(STRUCT,CLASS,PATH) \
template<> std::string tpcrs::ConfigStruct<St_ ## STRUCT ## C, CLASS , STRUCT ## _st>::name(# PATH);

#define MakeChairAltInstance(STRUCT,PATHA,PATHB,AorB)	\
template<> std::string tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_ ## STRUCT ## C, STRUCT ## _st>::name(AorB==0 ? # PATHA : # PATHB);

#define MakeChairAltInstance2(STRUCT,CLASS,PATHA,PATHB,AorB)	\
template<> std::string tpcrs::ConfigStruct<St_ ## STRUCT ## C, CLASS , STRUCT ## _st>::name(AorB==0 ? # PATHA : # PATHB);

static int _debug = 0;
//__________________Calibrations/tpc______________________________________________________________
#include "StDetectorDbMaker/St_tpcGasC.h"
MakeChairInstance(tpcGas,Calibrations/tpc/tpcGas);
#include "St_TpcEffectivedXC.h"
MakeChairInstance(TpcEffectivedX,Calibrations/tpc/TpcEffectivedX);
#include "StDetectorDbMaker/St_tpcGridLeakC.h"
MakeChairInstance(tpcGridLeak,Calibrations/tpc/tpcGridLeak);
#include "StDetectorDbMaker/St_tpcOmegaTauC.h"
MakeChairInstance(tpcOmegaTau,Calibrations/tpc/tpcOmegaTau);
#include "StDetectorDbMaker/St_tpcDriftVelocityC.h"
MakeChairInstance(tpcDriftVelocity,Calibrations/tpc/tpcDriftVelocity);
#include "St_TpcSecRowBC.h"
MakeChairInstance2(TpcSecRowCor,St_TpcSecRowBC,Calibrations/tpc/TpcSecRowB);
#include "St_TpcSecRowCC.h"
MakeChairInstance2(TpcSecRowCor,St_TpcSecRowCC,Calibrations/tpc/TpcSecRowC);
#include "StDetectorDbMaker/St_tpcCorrectionC.h"
#include "StDetectorDbMaker/St_tpcCalibResolutionsC.h"
MakeChairInstance(tpcCalibResolutions,Calibrations/tpc/tpcCalibResolutions);
#include "StDetectorDbMaker/St_tpcChargeEventC.h"
MakeChairInstance(tpcChargeEvent,Calibrations/tpc/tpcChargeEvent);
#include "StDetectorDbMaker/St_tpcSCGLC.h"
MakeChairInstance(tpcSCGL,Calibrations/tpc/tpcSCGL);


double St_tpcCorrectionC::CalcCorrection(int i, double x, double z, int NparMax) {
  tpcCorrection_st *cor =  Struct(i);
  return SumSeries(cor, x, z, NparMax);
}


double St_tpcCorrectionC::SumSeries(tpcCorrection_st *cor,  double x, double z, int NparMax) {
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
#include "St_TpcRowQC.h"
MakeChairInstance2(tpcCorrection,St_TpcRowQC,Calibrations/tpc/TpcRowQ);
#include "St_TpcDriftDistOxygenC.h"
MakeChairInstance2(tpcCorrection,St_TpcDriftDistOxygenC,Calibrations/tpc/TpcDriftDistOxygen);
#include "St_TpcMultiplicityC.h"
MakeChairInstance2(tpcCorrection,St_TpcMultiplicityC,Calibrations/tpc/TpcMultiplicity);
#include "St_TpcZCorrectionBC.h"
MakeChairInstance2(tpcCorrection,St_TpcZCorrectionBC,Calibrations/tpc/TpcZCorrectionB);
#include "St_TpcdXCorrectionBC.h"
MakeChairInstance2(tpcCorrection,St_TpcdXCorrectionBC,Calibrations/tpc/TpcdXCorrectionB);
#include "St_tpcPressureBC.h"
MakeChairInstance2(tpcCorrection,St_tpcPressureBC,Calibrations/tpc/tpcPressureB);
#include "St_TpcEdgeC.h"
MakeChairInstance2(tpcCorrection,St_TpcEdgeC,Calibrations/tpc/TpcEdge);
#include "St_TpcAdcCorrectionBC.h"
MakeChairInstance2(tpcCorrection,St_TpcAdcCorrectionBC,Calibrations/tpc/TpcAdcCorrectionB);
#include "St_TpcAdcCorrectionMDF.h"
MakeChairInstance2(MDFCorrection,St_TpcAdcCorrectionMDF,Calibrations/tpc/TpcAdcCorrectionMDF);
#include "St_tpcMethaneInC.h"
MakeChairInstance2(tpcCorrection,St_tpcMethaneInC,Calibrations/tpc/tpcMethaneIn);
#include "St_tpcGasTemperatureC.h"
MakeChairInstance2(tpcCorrection,St_tpcGasTemperatureC,Calibrations/tpc/tpcGasTemperature);
#include "St_tpcWaterOutC.h"
MakeChairInstance2(tpcCorrection,St_tpcWaterOutC,Calibrations/tpc/tpcWaterOut);
#include "St_tpcTimeDependenceC.h"
MakeChairInstance2(tpcCorrection,St_tpcTimeDependenceC,Calibrations/tpc/tpcTimeDependence);
#include "St_TpcdChargeC.h"
MakeChairInstance2(tpcCorrection,St_TpcdChargeC,Calibrations/tpc/TpcdCharge);
#include "St_TpcrChargeC.h"
MakeChairInstance2(tpcCorrection,St_TpcrChargeC,Calibrations/tpc/TpcrCharge);
#include "St_TpcTanLC.h"
MakeChairInstance2(tpcCorrection,St_TpcTanLC,Calibrations/tpc/TpcTanL);
#include "St_TpcCurrentCorrectionC.h"
MakeChairInstance2(tpcCorrection,St_TpcCurrentCorrectionC,Calibrations/tpc/TpcCurrentCorrectionX);
#include "St_TpcZDCC.h"
MakeChairInstance2(tpcCorrection,St_TpcZDCC,Calibrations/tpc/TpcZDC);
#include "St_TpcSpaceChargeC.h"
MakeChairInstance2(tpcCorrection,St_TpcSpaceChargeC,Calibrations/tpc/TpcSpaceCharge);
#include "St_TpcPhiDirectionC.h"
MakeChairInstance2(tpcCorrection,St_TpcPhiDirectionC,Calibrations/tpc/TpcPhiDirection);
#include "St_TpcdEdxCorC.h"
MakeChairInstance2(tpcCorrection,St_TpcdEdxCorC,Calibrations/tpc/TpcdEdxCor);
#include "St_TpcLengthCorrectionBC.h"
MakeChairInstance2(tpcCorrection,St_TpcLengthCorrectionBC,Calibrations/tpc/TpcLengthCorrectionB);
#include "St_TpcLengthCorrectionMDF.h"
MakeChairInstance2(MDFCorrection,St_TpcLengthCorrectionMDF,Calibrations/tpc/TpcLengthCorrectionMDF);
#include "St_TpcPadCorrectionMDF.h"
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
  double returnValue = fFunc[k]->GetSave(xx);
  return returnValue;
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
      double y  =  1 + 2. / (XMax(k)[j] - XMin(k)[j])
	* (x[j] - XMax(k)[j]);
      term        *= EvalFactor(p,y);
    }
    // Add this term to the final result
    returnValue += term*term;
  }
  returnValue = std::sqrt(returnValue);
  return returnValue;
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
#include "StDetectorDbMaker/St_tpcEffectiveGeomC.h"
MakeChairAltInstance(tpcEffectiveGeom,Calibrations/tpc/tpcEffectiveGeom,Calibrations/tpc/tpcEffectiveGeomB,gEnv->GetValue("NewTpcAlignment",0));
#include "StDetectorDbMaker/St_tpcElectronicsC.h"
MakeChairAltInstance(tpcElectronics,Calibrations/tpc/tpcElectronics,Calibrations/tpc/tpcElectronicsB,gEnv->GetValue("NewTpcAlignment",0));
#include "StDetectorDbMaker/St_tpcPadResponseC.h"
MakeChairInstance(tpcPadResponse,Calibrations/tpc/tpcPadResponse);
#include "StDetectorDbMaker/St_tpcHighVoltagesC.h"
MakeChairInstance(tpcHighVoltages,Calibrations/tpc/tpcHighVoltages);
#include "StDetectorDbMaker/St_tpcPadrowT0C.h"
MakeChairAltInstance(tpcPadrowT0,Calibrations/tpc/tpcPadrowT0,Calibrations/tpc/tpcPadrowT0B,gEnv->GetValue("NewTpcAlignment",0));
#include "StDetectorDbMaker/St_tpcSectorT0offsetC.h"
MakeChairInstance(tpcSectorT0offset,Calibrations/tpc/tpcSectorT0offset);
#include "StDetectorDbMaker/St_tpcAltroParamsC.h"
MakeChairInstance(tpcAltroParams,Calibrations/tpc/tpcAltroParams);
#include "StDetectorDbMaker/St_asic_thresholdsC.h"
MakeChairInstance(asic_thresholds,Calibrations/tpc/asic_thresholds);
#include "StDetectorDbMaker/St_tpcAnodeHVC.h"
MakeChairInstance(tpcAnodeHV,Calibrations/tpc/tpcAnodeHV);
#include "StDetectorDbMaker/St_tpcPadPlanesC.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
unsigned char          St_tpcPadConfigC::iTpc(int sector)                     {unsigned char iTPC = Struct()->itpc[sector-1];  return iTPC;}
int 	         St_tpcPadConfigC::padRows(int sector) 	          {return St_tpcPadPlanesC::instance()->padRows()               ;}
int 	         St_tpcPadConfigC::innerPadRows(int sector) 	          {return St_tpcPadPlanesC::instance()->innerPadRows() 	   ;}
int 	         St_tpcPadConfigC::innerPadRows48(int sector) 	  {return St_tpcPadPlanesC::instance()->innerPadRows48()	   ;}
int 	         St_tpcPadConfigC::innerPadRows52(int sector) 	  {return St_tpcPadPlanesC::instance()->innerPadRows52()	   ;}
int 	         St_tpcPadConfigC::outerPadRows(int sector) 	          {return St_tpcPadPlanesC::instance()->outerPadRows() 	   ;}
int 	         St_tpcPadConfigC::superInnerPadRows(int sector)        {return St_tpcPadPlanesC::instance()->superInnerPadRows()     ;}
int 	         St_tpcPadConfigC::superOuterPadRows(int sector)        {return St_tpcPadPlanesC::instance()->superOuterPadRows()     ;}
double 	 St_tpcPadConfigC::innerSectorPadWidth(int sector)      {return St_tpcPadPlanesC::instance()->innerSectorPadWidth()   ;}
double 	 St_tpcPadConfigC::innerSectorPadLength(int sector)     {return St_tpcPadPlanesC::instance()->innerSectorPadLength()  ;}
double 	 St_tpcPadConfigC::innerSectorPadPitch(int sector)      {return St_tpcPadPlanesC::instance()->innerSectorPadPitch()   ;}
double 	 St_tpcPadConfigC::innerSectorRowPitch1(int sector)     {return St_tpcPadPlanesC::instance()->innerSectorRowPitch1()  ;}
double 	 St_tpcPadConfigC::innerSectorRowPitch2(int sector)     {return St_tpcPadPlanesC::instance()->innerSectorRowPitch2()  ;}
double 	 St_tpcPadConfigC::firstPadRow(int sector) 	          {return St_tpcPadPlanesC::instance()->firstPadRow() 	   ;}
double 	 St_tpcPadConfigC::firstOuterSectorPadRow(int sector)   {return St_tpcPadPlanesC::instance()->firstOuterSectorPadRow();}
double 	 St_tpcPadConfigC::lastOuterSectorPadRow(int sector)    {return St_tpcPadPlanesC::instance()->lastOuterSectorPadRow() ;}
double 	 St_tpcPadConfigC::firstRowWidth(int sector) 	          {return St_tpcPadPlanesC::instance()->firstRowWidth()         ;}
double 	 St_tpcPadConfigC::lastRowWidth(int sector) 	          {return St_tpcPadPlanesC::instance()->lastRowWidth()          ;}
double 	 St_tpcPadConfigC::outerSectorPadWidth(int sector)      {return St_tpcPadPlanesC::instance()->outerSectorPadWidth()   ;}
double 	 St_tpcPadConfigC::outerSectorPadLength(int sector)     {return St_tpcPadPlanesC::instance()->outerSectorPadLength()  ;}
double 	 St_tpcPadConfigC::outerSectorPadPitch(int sector)      {return St_tpcPadPlanesC::instance()->outerSectorPadPitch()   ;}
double 	 St_tpcPadConfigC::outerSectorRowPitch(int sector)      {return St_tpcPadPlanesC::instance()->outerSectorRowPitch()   ;}
double 	 St_tpcPadConfigC::outerSectorLength(int sector)        {return St_tpcPadPlanesC::instance()->outerSectorLength()     ;}
double 	 St_tpcPadConfigC::ioSectorSeparation(int sector)       {return St_tpcPadPlanesC::instance()->ioSectorSeparation()    ;}
double 	 St_tpcPadConfigC::innerSectorEdge(int sector) 	  {return St_tpcPadPlanesC::instance()->innerSectorEdge()       ;}
double 	 St_tpcPadConfigC::outerSectorEdge(int sector) 	  {return St_tpcPadPlanesC::instance()->outerSectorEdge() 	   ;}
double 	 St_tpcPadConfigC::innerSectorPadPlaneZ(int sector)     {return St_tpcPadPlanesC::instance()->innerSectorPadPlaneZ()  ;}
double 	 St_tpcPadConfigC::outerSectorPadPlaneZ(int sector)     {return St_tpcPadPlanesC::instance()->outerSectorPadPlaneZ()  ;}
int* 	         St_tpcPadConfigC::innerPadsPerRow(int sector) 	  {return St_tpcPadPlanesC::instance()->innerPadsPerRow()       ;}
int* 	         St_tpcPadConfigC::outerPadsPerRow(int sector) 	  {return St_tpcPadPlanesC::instance()->outerPadsPerRow()       ;}
int            St_tpcPadConfigC::padsPerRow(int sector, int row)    {
  int Ninner = innerPadRows(sector);
  return (row <= Ninner) ?
    innerPadsPerRow(sector)[row-1] :
    outerPadsPerRow(sector)[row-1-Ninner];
}
double* 	 St_tpcPadConfigC::innerRowRadii(int sector) 	          {return St_tpcPadPlanesC::instance()->innerRowRadii()         ;}
double* 	 St_tpcPadConfigC::outerRowRadii(int sector) 	          {return St_tpcPadPlanesC::instance()->outerRowRadii() 	   ;}
// taken from StRItpcPadPlane
int            St_tpcPadConfigC::numberOfRows(int sector)             {return padRows(sector);}
int            St_tpcPadConfigC::numberOfInnerRows(int sector)        {return innerPadRows(sector);}
int            St_tpcPadConfigC::numberOfInnerRows48(int sector)      {return innerPadRows48(sector);}
int            St_tpcPadConfigC::numberOfInnerRows52(int sector)      {return innerPadRows52(sector);}
int            St_tpcPadConfigC::numberOfOuterRows(int sector)        {return outerPadRows(sector);}
bool           St_tpcPadConfigC::isRowInRange(int sector, int row)  {return (row >= 1 && row<=numberOfRows(sector)) ? true: false;}
double         St_tpcPadConfigC::radialDistanceAtRow(int sector, int row)       {
  if (! isRowInRange(sector,row)) return 0;
  int Ninner = innerPadRows(sector);
  if ( row<=Ninner ) return innerRowRadii(sector)[row-1];
  else               return outerRowRadii(sector)[row-1-Ninner];
}
int            St_tpcPadConfigC::numberOfPadsAtRow(int sector, int row)       {
  if (! isRowInRange(sector, row)) return 0;
  int Ninner = innerPadRows(sector);
  if ( row<=Ninner ) return innerPadsPerRow(sector)[row-1];
  return outerPadsPerRow(sector)[row-1-Ninner];
}
double         St_tpcPadConfigC::PadWidthAtRow(int sector, int row)       {
  if (! isRowInRange(sector,row)) return 0;
  int Ninner = innerPadRows(sector);
  if ( row<=Ninner) return innerSectorPadWidth(sector);
  return outerSectorPadWidth(sector);
}
double         St_tpcPadConfigC::PadLengthAtRow(int sector, int row)       {
  if (! isRowInRange(sector,row)) return 0;
  int Ninner = innerPadRows(sector);
  if ( row<=Ninner) return innerSectorPadLength(sector);
  return outerSectorPadLength(sector);
}
double         St_tpcPadConfigC::PadPitchAtRow(int sector, int row)       {
  if (! isRowInRange(sector,row)) return 0;
  int Ninner = innerPadRows(sector);
  if ( row<=Ninner) return innerSectorPadPitch(sector);
  return outerSectorPadPitch(sector);
}
double         St_tpcPadConfigC::RowPitchAtRow(int sector, int row)       {
  if (! isRowInRange(sector,row)) return 0;
  int Ninner = innerPadRows(sector);
  if ( row<=numberOfInnerRows48(sector) ) return innerSectorRowPitch1(sector);
  else if (row>numberOfInnerRows48(sector)&&row<=Ninner) return innerSectorRowPitch2(sector);
  return outerSectorRowPitch(sector);
}
int            St_tpcPadConfigC::indexForRowPad(int sector, int row, int pad)       {
  if (pad >numberOfPadsAtRow(sector,row)) return -1;
  int index = 0;
  int Ninner = innerPadRows(sector);
  if (row>0 && row<=Ninner )             for (int i=1;i<row;i++) index += numberOfPadsAtRow(sector,i);
  else
    if (row>Ninner&&row<=numberOfRows(sector)) for (int i=Ninner+1;i<row;i++)  index += numberOfPadsAtRow(sector,i);
  index+=pad-1;
  return index;
}
#include "StDetectorDbMaker/St_TpcAvgPowerSupplyC.h"


void  St_tpcAnodeHVC::sockets(int sector, int padrow, int &e1, int &e2, float &f2) {
  e1 = (sector-1)*19;
  e2 = e1;
  f2 = 0;
  // sector=1..24 , padrow=1..45
  // f2 represents signal couplings from neighboring HV sections
  // see: http://www.star.bnl.gov/public/tpc/hard/signals/signal_division.html
  if (!  St_tpcPadConfigC::instance()->iTPC(sector)) {
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
  } else { // iTPC
    switch (padrow) {
    case  1:
    case  2:
    case  3: e1+= 1; e2+= 1; break;
    case  4: e1+= 1; e2+= 2; break;
    case  5:
    case  6:
    case  7:
    case  8: e1+= 2; e2+= 2; break;
    case  9: e1+= 2; e2+= 3; break;
    case 10: 
    case 11:
    case 12:
    case 13: e1+= 3; e2+= 3; break;
    case 14: e1+= 3; e2+= 4; break;
    case 15:
    case 16:
    case 17:
    case 18: e1+= 4; e2+= 4; break;
    case 19: e1+= 4; e2+= 5; break;
    case 20:
    case 21:
    case 22:
    case 23: e1+= 5; e2+= 5; break;
    case 24: e1+= 5; e2+= 6; break;
    case 25:
    case 26:
    case 27:
    case 28: e1+= 6; e2+= 6; break;
    case 29: e1+= 6; e2+= 7; break;
    case 30:
    case 31:
    case 32:
    case 33: e1+= 7; e2+= 7; break;
    case 34: e1+= 7; e2+= 8; break;
    case 35:
    case 36:
    case 37:
    case 38: e1+= 8; e2+= 8; break;
    case 39:
    case 40: e1+=17; e2+=17; break;
    case 41: e1+= 9; e2+=18;  break;
    case 42:
    case 43: e1+= 9; break;
    case 44: e1+= 9; e2+=10;  break;
    case 45:
    case 46:
    case 47: e1+=10; break;
    case 48: e1+=10; e2+=11;  break;
    case 49:
    case 50:
    case 51: e1+=11; break;
    case 52: e1+=11; e2+=12; f2 = 0.40250; break;
    case 53:
    case 54:
    case 55: e1+=12; break;
    case 56: e1+=12; e2+=13; f2 = 0.40250; break;
    case 57:
    case 58:
    case 59: e1+=13; break;
    case 60: e1+=13; e2+=14; f2 = 0.40250; break;
    case 61:
    case 62:
    case 63: e1+=14; break;
    case 64: e1+=14; e2+=15; f2 = 0.40250; break;
    case 65:
    case 66:
    case 67: e1+=15; break;
    case 68: e1+=15; e2+=16; f2 = 0.40250; break;
    case 69:
    case 70:
    case 71: e1+=16; break;
    case 72: e1+=16; e2+=19; f2 = 0.40250; break;
    default: e1 = 0; e2 = 0; f2 = 0;
    }
  }
}


float St_tpcAnodeHVC::voltage(int i) const {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
     LOG_ERROR << "St_tpcAnodeHVC::voltage(" << i << " is called but the valid St_TpcAvgPowerSupplyC::instance() exists\n";
  }
  return Struct(i)->voltage;
}


float St_tpcAnodeHVC::voltagePadrow(int sector, int padrow) const {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->voltagePadrow(sector,padrow);
  }
  int e1 = 0, e2 = 0;
  float f2 = 0;
  St_tpcAnodeHVC::sockets(sector, padrow, e1, e2, f2);
  if (e1==0) return -99;
  float v1=voltage(e1-1);
  if (f2 < 0.1) return v1;
  float v2=voltage(e2-1);
  if (std::abs(v2 - v1) > 40) return -99;
  if (std::abs(v2 - v1) <  1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  float B = (padrow <= St_tpcPadConfigC::instance()->innerPadRows(sector) ? 13.05e-3 : 10.26e-3);
  float v_eff = std::log((1.0-f2)*std::exp(B*v1) + f2*std::exp(B*v2)) / B;
  return v_eff;
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
  float B = (padrow <= St_tpcPadConfigC::instance()->innerPadRows(sector) ? 13.05e-3 : 10.26e-3);
  float v_eff = std::log((1.0-f2)*std::exp(B*v1) + f2*std::exp(B*v2)) / B;
  return v_eff;
}


float St_TpcAvgPowerSupplyC::AcChargeL(int sector, int channel) {
  //  static const double RA[2]        = { 154.484, 81.42}; // Outer/ Inner average Radii
  //  static const double WireLenth[2] = {   3.6e5, 1.6e5};
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

#include "StDetectorDbMaker/St_tpcAnodeHVavgC.h"
MakeChairInstance(tpcAnodeHVavg,Calibrations/tpc/tpcAnodeHVavg);


float St_tpcAnodeHVavgC::voltage(int i) const {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
     LOG_ERROR << "St_tpcAnodeHVavgC::voltage(" << i << " is called but the valid St_TpcAvgPowerSupplyC::instance() exists\n";
  }
  return Struct(i)->voltage;
}


bool St_tpcAnodeHVavgC::tripped(int sector, int padrow) const {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->tripped(sector,padrow);
  }
  return (voltage() < -100);
}


float St_tpcAnodeHVavgC::voltagePadrow(int sector, int padrow) const {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->voltagePadrow(sector,padrow);
  }
  int e1 = 0, e2 = 0;
  float f2 = 0;
  St_tpcAnodeHVC::sockets(sector, padrow, e1, e2, f2);
  if (e1==0) return -99;
  float v1=voltage(e1-1);
  if (f2==0) return v1;
  float v2=voltage(e2-1);
  if (v2==v1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  float B = (padrow <= St_tpcPadConfigC::instance()->innerPadRows(sector) ? 13.05e-3 : 10.26e-3);
  float v_eff = std::log((1.0-f2)*std::exp(B*v1) + f2*std::exp(B*v2)) / B;
  return v_eff;
}


#include "StDetectorDbMaker/St_tpcPadGainT0C.h"
MakeChairInstance(tpcPadGainT0,Calibrations/tpc/tpcPadGainT0);
#include "StDetectorDbMaker/St_itpcPadGainT0C.h"
MakeChairInstance(itpcPadGainT0,Calibrations/tpc/itpcPadGainT0);
#include "StDetectorDbMaker/St_tpcPadGainT0BC.h"
St_tpcPadGainT0BC *St_tpcPadGainT0BC::fgInstance = 0;
St_tpcPadGainT0BC *St_tpcPadGainT0BC::instance() {if (! fgInstance) fgInstance = new St_tpcPadGainT0BC(); return fgInstance;}


float 	St_tpcPadGainT0BC::Gain(int sector, int row, int pad) const {
  float gain = 0;
  if (St_tpcPadConfigC::instance()->iTPC(sector)) {
    if (row <= 40) {
      gain = St_itpcPadGainT0C::instance()->Gain(sector,row,pad);
    } else {
      gain = St_tpcPadGainT0C::instance()->Gain(sector,row-40+13,pad);
    }
  } else { // Tpx
    gain = St_tpcPadGainT0C::instance()->Gain(sector,row,pad);
  }
  return gain;
}


float 	  St_tpcPadGainT0BC::T0(int sector, int row, int pad) const {
  float T0 = 0;
  if (St_tpcPadConfigC::instance()->iTPC(sector)) {
    if (row <= 40) 
      T0 = St_itpcPadGainT0C::instance()->T0(sector,row,pad);
    else 
      T0 = St_tpcPadGainT0C::instance()->T0(sector,row-40+13,pad);
  } else { // Tpx
    T0 = St_tpcPadGainT0C::instance()->T0(sector,row,pad);
  }
  return T0;
}


bool    St_tpcPadGainT0BC::livePadrow(int sector, int row) const {
  if (St_tpcPadConfigC::instance()->iTPC(sector)) {
    if (row <= 40)
      return St_itpcPadGainT0C::instance()->livePadrow(sector,row);
    else 
      return St_tpcPadGainT0C::instance()->livePadrow(sector,row-40+13);
  }
  return St_tpcPadGainT0C::instance()->livePadrow(sector,row);
}


#include "StDetectorDbMaker/St_TpcResponseSimulatorC.h"
MakeChairInstance(TpcResponseSimulator,Calibrations/tpc/TpcResponseSimulator);
#include "StDetectorDbMaker/St_tpcGainCorrectionC.h"
MakeChairInstance2(tpcCorrection,St_tpcGainCorrectionC,Calibrations/tpc/tpcGainCorrection);
#include "StDetectorDbMaker/St_TpcAvgCurrentC.h"
MakeChairInstance(TpcAvgCurrent,Calibrations/tpc/TpcAvgCurrent);


int St_TpcAvgCurrentC::ChannelFromRow(int sector, int row) {
  if (row <  1 || row > St_tpcPadConfigC::instance()->padRows(sector)) return -1;
  if (!  St_tpcPadConfigC::instance()->iTPC(sector)) {
    if (row <  3) return 1;
    if (row <  7) return 2;
    if (row < 10) return 3;
    if (row < 14) return 4;
    if (row < 22) return 5;
    if (row < 30) return 6;
    if (row < 38) return 7;
    return 8;
  } else { // iTPC
    // Jim Thomas, mail from 09/27/17
    if (row < 10) return 1; //  9 shared 1&2
    if (row < 20) return 2; // 19 shared 2&3
    if (row < 30) return 3; // 29 shared 3&4
    if (row < 14 - 13 + 40) return 4;
    if (row < 22 - 13 + 40) return 5;
    if (row < 30 - 13 + 40) return 6;
    if (row < 38 - 13 + 40) return 7;
    return 8;
  }
  return -1;
}


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
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->AcChargeL(sector,channel);
  }
  //  static const double RA[2]        = { 154.484, 81.42}; // Outer/ Inner average Radii
  //  static const double WireLenth[2] = {   3.6e5, 1.6e5};
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


float St_TpcAvgCurrentC::AvCurrent(int sector, int channel) {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->AvCurrent(sector,channel);
  }
  return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
    Struct()->AvCurrent[8*(sector-1)+channel-1] :     0;
}


float St_TpcAvgCurrentC::AcCharge(int sector, int channel) {
  if (! St_TpcAvgPowerSupplyC::instance()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->AcCharge(sector,channel);
  }
  return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
    Struct()->AcCharge[8*(sector-1)+channel-1] :     0;
}
#include "StDetectorDbMaker/St_tpcRDOMapC.h"
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
#include "StDetectorDbMaker/St_tpcRDOT0offsetC.h"
MakeChairInstance(tpcRDOT0offset,Calibrations/tpc/tpcRDOT0offset);
float St_tpcRDOT0offsetC::T0(int sector, int padrow, int pad) const {
  float t0 = 0;
  if (! IsShfited(sector)) return t0;
  if (St_tpcPadConfigC::instance()->iTPC(sector) && padrow <= 40)  return t0; // no shift in iTPC
  int rdo = St_tpcRDOMapC::instance()->rdo(padrow,pad);
  if (!rdo) return t0;
  t0 = Struct()->t0[sector-1][rdo-1];
  return t0;
}

#include "StDetectorDbMaker/St_trigDetSumsC.h"
MakeChairInstance(trigDetSums, Calibrations/rich/trigDetSums);
//___________________tpc_____________________________________________________________
#include "StDetectorDbMaker/St_tss_tssparC.h"
MakeChairInstance(tss_tsspar,tpc/tsspars/tsspar);


float St_tss_tssparC::gain(int sector, int row) {
  int l = 0;
  double V_nominal = 1390;
  float V = 0;
  float gain = 0;
  if (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) {l = 1; V_nominal = 1170;}
  St_tpcGainCorrectionC *gC = St_tpcGainCorrectionC::instance();
  int NRows = gC->GetNRows();
  if (l >= NRows) return gain;
  V = St_tpcAnodeHVavgC::instance()->voltagePadrow(sector,row);
  if (V > 0) {
    double v = V - V_nominal;
    if (v < gC->min(l) || v > gC->max(l)) return gain;
    if (gC->min(l) < -450) {
      // if range was expanded below 150 V then use only the linear approximation
      gain  = std::exp(gC->CalcCorrection(l,v, 0., 2));
    } else {
      gain  = std::exp(gC->CalcCorrection(l,v));
    }
  }
  return gain;
}
//__________________Calibrations/rich______________________________________________________________
#include "StDetectorDbMaker/St_richvoltagesC.h"
MakeChairInstance(richvoltages,Calibrations/rich/richvoltages);
#include "StDetectorDbMaker/St_spaceChargeCorC.h"
MakeChairInstance2(spaceChargeCor,St_spaceChargeCorR1C,Calibrations/rich/spaceChargeCor);
MakeChairInstance2(spaceChargeCor,St_spaceChargeCorR2C,Calibrations/rich/spaceChargeCorR2);
//_________________RunLog_____________________________________________________________
#include "St_MagFactorC.h"
MakeChairInstance(MagFactor,RunLog/MagFactor);
//_________________RunLog/onl_______________________________________________________________
#include "StDetectorDbMaker/St_starClockOnlC.h"
MakeChairInstance(starClockOnl,RunLog/onl/starClockOnl);

#include "StDetectorDbMaker/St_starMagOnlC.h"
MakeChairInstance(starMagOnl,RunLog/onl/starMagOnl);

#include "StDetectorDbMaker/St_tpcRDOMasksC.h"
MakeChairInstance(tpcRDOMasks,RunLog/onl/tpcRDOMasks);


unsigned int       St_tpcRDOMasksC::getSectorMask(unsigned int sector) {
  unsigned int MASK = 0x0000; // default is to mask it out
  //unsigned int MASK = 0xFFFF; // change to  ON by default ** THIS WAS A HACK
  if(sector < 1 || sector > 24 || GetNRows() == 0){
    LOG_WARN << "St_tpcRDOMasksC:: getSectorMask : return default mask for "
     << "sector= " << sector << " GetNRows()=" << GetNRows() << '\n';
    return MASK;
  }
  MASK = mask(((sector + 1) / 2) - 1); // does the mapping from sector 1-24 to packed sectors
  if( sector % 2 == 0){ // if its even relevent bits are 6-11
    MASK = MASK >> 6;
  }
  // Otherwise want lower 6 bits
  MASK &= 0x000003F; // Mask out higher order bits
  if (sector == 16 && MASK == 0 && runNumber() > 8181000 && runNumber() < 9181000) MASK = 4095;
  return MASK;
}

//___________________Conditions/trg_____________________________________________________________
#include "StDetectorDbMaker/St_trgTimeOffsetC.h"
MakeChairAltInstance(trgTimeOffset,Conditions/trg/trgTimeOffset,Conditions/trg/trgTimeOffsetB,gEnv->GetValue("NewTpcAlignment",0));
//___________________Geometry/tpc_____________________________________________________________
#include "StDetectorDbMaker/St_tpcDimensionsC.h"
MakeChairInstance(tpcDimensions,Geometry/tpc/tpcDimensions);
#include "StDetectorDbMaker/St_tpcWirePlanesC.h"
MakeChairInstance(tpcWirePlanes,Geometry/tpc/tpcWirePlanes);
#include "StDetectorDbMaker/St_tpcFieldCageC.h"
MakeChairInstance(tpcFieldCage,Geometry/tpc/tpcFieldCage);
MakeChairInstance(tpcPadPlanes,Geometry/tpc/tpcPadPlanes);
MakeChairInstance(tpcPadConfig,Geometry/tpc/tpcPadConfig);
#include "StDetectorDbMaker/St_tpcGlobalPositionC.h"
MakeChairInstance(tpcGlobalPosition,Geometry/tpc/tpcGlobalPosition);
#include "StDetectorDbMaker/St_tpcFieldCageShortC.h"
MakeChairInstance(tpcFieldCageShort,Geometry/tpc/tpcFieldCageShort);
#include "StDetectorDbMaker/St_tpcHVPlanesC.h"
MakeChairInstance(tpcHVPlanes,Geometry/tpc/tpcHVPlanes);
#include "StDetectorDbMaker/St_SurveyC.h"
#include "StDetectorDbMaker/StTpcSurveyC.h"
MakeChairAltInstance2(Survey,StTpcInnerSectorPosition,Geometry/tpc/TpcInnerSectorPosition,Geometry/tpc/TpcInnerSectorPositionB,gEnv->GetValue("NewTpcAlignment",0));
MakeChairAltInstance2(Survey,StTpcOuterSectorPosition,Geometry/tpc/TpcOuterSectorPosition,Geometry/tpc/TpcOuterSectorPositionB,gEnv->GetValue("NewTpcAlignment",0));
MakeChairAltInstance2(Survey,StTpcSuperSectorPosition,Geometry/tpc/TpcSuperSectorPosition,Geometry/tpc/TpcSuperSectorPositionB,gEnv->GetValue("NewTpcAlignment",0));
MakeChairInstance2(Survey,StTpcHalfPosition,Geometry/tpc/TpcHalfPosition);
MakeChairInstance2(Survey,StTpcPosition,Geometry/tpc/TpcPosition);
#include "StDetectorDbMaker/St_iTPCSurveyC.h"
MakeChairInstance(iTPCSurvey,Geometry/tpc/iTPCSurvey);

St_SurveyC::St_SurveyC() : fRotations(0)  { }


St_SurveyC::~St_SurveyC() {
  if (fRotations) {
    for (unsigned int i = 0; i < GetNRows(); i++) {
      SafeDelete(fRotations[0]);
    }
    SafeDelete(fRotations);
  }
}


double St_SurveyC::IsOrtogonal(const double *r) {
// Perform orthogonality test for rotation.
  double cmax = 0;
  double cij;
  for (int i=0; i<2; i++) {
    for (int j=i+1; j<3; j++) {
      // check columns
      cij = std::abs(r[i]*r[j]+r[i+3]*r[j+3]+r[i+6]*r[j+6]);
      if (cij>1E-4) cmax = cij;
      // check rows
      cij = std::abs(r[3*i]*r[3*j]+r[3*i+1]*r[3*j+1]+r[3*i+2]*r[3*j+2]);
      if (cij>cmax) cmax = cij;
    }
  }
  return cmax;
}


void St_SurveyC::Normalize(TGeoHMatrix &R) {
  if (_debug) {
    LOG_INFO << "Matrix:\n"; R.Print("");
    LOG_INFO << "Determinant-1 = " << R.Determinant()-1 << '\n';
    const double *rr = R.GetRotationMatrix();
    LOG_INFO << "Ortogonality " << IsOrtogonal(rr) << '\n';
  }
  return;
}


const TGeoHMatrix &St_SurveyC::GetMatrix(int i) {
  assert(fRotations || fRotations[i]);
  assert(std::abs(fRotations[i]->Determinant())-1 < 1.e-3);
  return *fRotations[i];
}


const TGeoHMatrix &St_SurveyC::GetMatrix4Id(int id) {
  for (unsigned int i = 0; i < GetNRows(); i++) {
    if (Id(i) == id) {
      return GetMatrix(i);
    }
  }
  LOG_INFO  << "St_SurveyC::GetMatrix4Id(" << id << ") entry has not been found\n";
  assert(0);
  return GetMatrix(0);
}


const TGeoHMatrix &St_SurveyC::GetMatrixR(int i) {
  static TGeoHMatrix rot;
  double rotations[9] = {
    r00(i), r01(i),      0,
    r10(i), r11(i),      0,
    0     ,      0, r22(i)};
  rot.SetName(Form("%s_%i",GetName().c_str(),i));
  rot.SetRotation(rotations);
  rot.SetTranslation(Translation(i));
  return *&rot;
}


void St_SurveyC::GetAngles(double &phi, double &the, double &psi, int i) {
  phi = the = psi = 0;  // Korn 14.10-5
  double cosDelta = (r00(i) + r11(i) + r22(i) - 1)/2; // (Tr(R) - 1)/2
  double Delta = std::acos(cosDelta);
  if (Delta < 0) Delta += 2*M_PI;
  double sinDelta2 = std::sin(Delta/2);
  if (std::abs(sinDelta2) < 1.e-7) return;
  double c[3] = {
    (r21(i) - r12(i))/(2*sinDelta2), // a32-a23
    (r02(i) - r20(i))/(2*sinDelta2), // a13-a31
    (r10(i) - r01(i))/(2*sinDelta2)  // a21-a12
  };
  double u = std::atan2(c[0],c[1]);
  double v = std::atan(c[2]*std::tan(Delta/2));
  phi = (v - u)/2 - M_PI_2;
  psi = (v + u)/2 - M_PI_2;
  the = 2*std::atan2(c[0]*std::sin(v),c[2]*std::sin(u));
  double raddeg = 180./M_PI;
  phi   *= raddeg;
  the   *= raddeg;
  psi   *= raddeg;
}


double St_spaceChargeCorC::getSpaceChargeCoulombs(double scaleFactor)
  {
    St_trigDetSumsC* scalers = St_trigDetSumsC::instance();
    if (! scalers ) return 0;
    double zf = zeroField(0); // potential validity margin for scalers
    if (zf>0 && zf<1) scalers->setValidityMargin(zf);
    double coulombs = 0;

    bool use_powers = true;

    for (int row=0;row< (int) GetNRows();row++) {
      double mult = 0;
      switch ((int) getSpaceChargeDetector(row)) {
        case (0) : mult = scalers->getMult(); break; // vpdx as of 2007-12-19
        case (1) : mult = scalers->getBBCX(); break;
        case (2) : mult = scalers->getZDCX(); break;
        case (3) : mult = scalers->getZDCEast()+scalers->getZDCWest(); break;
        case (4) : mult = scalers->getBBCEast()+scalers->getBBCWest(); break;
        case (5) : mult = scalers->getZDCEast(); break;
        case (6) : mult = scalers->getZDCWest(); break;
        case (7) : mult = scalers->getBBCEast(); break;
        case (8) : mult = scalers->getBBCWest(); break;
        case (9) : mult = scalers->getBBCYellowBkg(); break;
        case (10): mult = scalers->getBBCBlueBkg(); break;
        case (11): mult = scalers->getPVPDEast(); break;
        case (12): mult = scalers->getPVPDWest(); break;
        case (13) : mult = scalers->getCTBOrTOFp(); break; // zdcx-no-killer as of 2011
        case (14) : mult = scalers->getCTBEast(); break; // zdce-no-killer as of 2011
        case (15) : mult = scalers->getCTBWest(); break; // zdcw-no-killer as of 2011

        default  : mult = 0.;
      }
      if (mult < 0) {
        Mark();
        return 0; // Unphysical scaler rates will be uncorrected
      } else UnMark();
      double saturation = getSpaceChargeSatRate(row);
      double correction = getSpaceChargeCorrection(scaleFactor,row);
      double factor     = getSpaceChargeFactor(row);
      double offset     = getSpaceChargeOffset(row);
      double intens = (mult < saturation) ? mult : saturation;
      if (use_powers) coulombs += ::pow(intens-offset,factor) * correction ;
      else coulombs += factor * (intens-offset) * correction ;
    }
    return coulombs;
  }

double St_tpcChargeEventC::timeDifference(unsigned long long bunchCrossingNumber, int idx) {
    // time difference between the event # idx and bunchCrossingNumber
    return ((double) (bunchCrossingNumber - eventBunchCrossing(idx))) /
      (St_starClockOnlC::instance()->Frequency());
  }

int St_tpcChargeEventC::indexBeforeBunchCrossing(unsigned long long bunchCrossingNumber) {
  // optimized search for typically looking for an index which
  // is the same or very close to previously found index

  int lastIndex = nChargeEvents() - 1;
  if (lastIndex < 0) return -999;

  if (localSearchUpperIndex < 0) { // no search yet
    // start from the middle and move outwards
    localSearchLowerIndex = nChargeEvents() / 2;
    localSearchUpperIndex = localSearchLowerIndex;
  }

  // try the same one as before first
  int direction = 0;
  if (bunchCrossingNumber >= eventBunchCrossing(localSearchUpperIndex) &&
      localSearchLowerIndex != lastIndex) direction = 1;
  else if (bunchCrossingNumber < eventBunchCrossing(localSearchLowerIndex)) direction = -1;
  else return localSearchLowerIndex;

  int delta = 1; // look up or down just one first

  while (direction > 0) { // look for higher indices
    localSearchLowerIndex = localSearchUpperIndex;
    localSearchUpperIndex = std::min(lastIndex,localSearchUpperIndex + delta);
    delta = localSearchUpperIndex - localSearchLowerIndex;
    if (bunchCrossingNumber < eventBunchCrossing(localSearchUpperIndex)) direction = 0;
    else if (localSearchUpperIndex == lastIndex) {
      localSearchLowerIndex = lastIndex;
      return lastIndex; // local indices will be lastIndex,lastIndex
    } else delta *= 2; // expand range and keep looking for higher indices
  }
  while (direction < 0) { // look for lower indices
    localSearchUpperIndex = localSearchLowerIndex;
    localSearchLowerIndex = std::max(0,localSearchLowerIndex - delta);
    delta = localSearchUpperIndex - localSearchLowerIndex;
    if (bunchCrossingNumber >= eventBunchCrossing(localSearchLowerIndex)) direction = 0;
    else if (localSearchLowerIndex == 0) {
      localSearchUpperIndex = 0;
      return -1; // local indices will be 0,0
    } else delta *= 2; // expand range and keep looking for lower indices
  }

  // already know that the result is within range
  while (delta > 1) {
    int tempIndex = localSearchLowerIndex + (delta/2);
    if (bunchCrossingNumber < eventBunchCrossing(tempIndex)) localSearchUpperIndex = tempIndex;
    else localSearchLowerIndex = tempIndex;
    delta = localSearchUpperIndex - localSearchLowerIndex;
  }
  // found
  return localSearchLowerIndex;
}

int St_tpcChargeEventC::findChargeTimes(unsigned long long bunchCrossingNumber, unsigned long long bunchCrossingWindow) {
  int idx2 = indexBeforeBunchCrossing(bunchCrossingNumber);
  int idx1 = indexBeforeBunchCrossing(bunchCrossingNumber-bunchCrossingWindow);
  int n = idx2-idx1;
  idx1++; // start from after the bunchCrossingWindow starts
  localStoreCharges.Set(n,&(eventCharges()[idx1]));
  localStoreTimesSinceCharges.Set(n);
  for (int i=0; i<n; i++) // must convert to times
    localStoreTimesSinceCharges.AddAt(timeDifference(bunchCrossingNumber,idx1+i),i);
  return n;
}

int St_tpcChargeEventC::findChargeTimes(unsigned long long bunchCrossingNumber, double timeWindow) {
  return findChargeTimes(bunchCrossingNumber,
    (unsigned long long) (timeWindow*(St_starClockOnlC::instance()->Frequency())));
}
