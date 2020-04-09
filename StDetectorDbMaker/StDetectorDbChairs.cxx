#include <cassert>
#include <string.h>
#include "TEnv.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "Math/SMatrix.h"
#include "tpcrs/logger.h"

using namespace ROOT::Math;

#include "St_base/StarChairDefs.h"

static Int_t _debug = 0;
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
MakeChairOptionalInstance2(TpcSecRowCor,St_TpcSecRowCC,Calibrations/tpc/TpcSecRowC);
#include "StDetectorDbMaker/St_tpcCorrectionC.h"
#include "StDetectorDbMaker/St_tpcCalibResolutionsC.h"
MakeChairInstance(tpcCalibResolutions,Calibrations/tpc/tpcCalibResolutions);
#include "StDetectorDbMaker/St_tpcChargeEventC.h"
MakeChairInstance(tpcChargeEvent,Calibrations/tpc/tpcChargeEvent);
#include "StDetectorDbMaker/St_tpcSCGLC.h"
MakeChairInstance(tpcSCGL,Calibrations/tpc/tpcSCGL);


Double_t St_tpcCorrectionC::CalcCorrection(Int_t i, Double_t x, Double_t z, Int_t NparMax) {
  tpcCorrection_st *cor =  ((St_tpcCorrection *) Table())->GetTable() + i;
  return SumSeries(cor, x, z, NparMax);
}


Double_t St_tpcCorrectionC::SumSeries(tpcCorrection_st *cor,  Double_t x, Double_t z, Int_t NparMax) {
  Double_t Sum = 0;
  if (! cor) return Sum;
  Int_t N = TMath::Abs(cor->npar)%100;
  if (N == 0) return Sum;
  if (NparMax > 0) N = NparMax;
  static Double_t T0, T1, T2;
  // parameterization variable
  Double_t X = x;
  if (cor->npar  < 0) X = TMath::Exp(x);
  else {
    switch  (cor->type) {
    case 10:// ADC correction offset + poly for ADC
    case 11:// ADC correction offset + poly for log(ADC) and |Z|
    case 12:// ADC correction offset + poly for log(ADC) and TanL
      X = TMath::Log(x);      break;
    case 1: // Tchebyshev [-1,1]
      if (cor->min < cor->max)   X = -1 + 2*TMath::Max(0.,TMath::Min(1.,(X - cor->min)/( cor->max - cor->min)));
      break;
    case 2: // Shifted TChebyshev [0,1]
      if (cor->min < cor->max)   X = TMath::Max(0.,TMath::Min(1.,(X - cor->min)/( cor->max - cor->min)));
      break;
    case 3:
      if (TMath::Abs(x) >= 1) X = 0;
      else                    X = TMath::Log(1. - TMath::Abs(x));
      break;
    case 4:
      if (TMath::Abs(x) >= 1) X = 0;
      else                    X = TMath::Sign(TMath::Log(1. - TMath::Abs(x)),x);
      break;
    case 5:
      if (x < 1e-7) X = -16.118;
      else          X = TMath::Log(x);
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
    Sum += TMath::Log(1. + cor->OffSet/x);
    Sum  = TMath::Exp(Sum);
    Sum *= x;
    break;
  case 11: // ADC correction offset + poly for log(ADC) and |Z|
    Sum = cor->a[1] + z*cor->a[2] + z*X*cor->a[3] + TMath::Exp(X*(cor->a[4] + X*cor->a[5]) + cor->a[6]);
    Sum *= TMath::Exp(-cor->a[0]);
    break;
  case 12: // ADC correction offset + poly for log(ADC) and TanL
    Sum = cor->a[1] + z*cor->a[2] + z*z*cor->a[3] + TMath::Exp(X*(cor->a[4] + X*cor->a[5]) + cor->a[6]);
    Sum *= TMath::Exp(-cor->a[0]);
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
MakeChairOptionalInstance2(tpcCorrection,St_TpcdChargeC,Calibrations/tpc/TpcdCharge);
#include "St_TpcrChargeC.h"
MakeChairOptionalInstance2(tpcCorrection,St_TpcrChargeC,Calibrations/tpc/TpcrCharge);
#include "St_TpcTanLC.h"
MakeChairInstance2(tpcCorrection,St_TpcTanLC,Calibrations/tpc/TpcTanL);
#include "St_TpcCurrentCorrectionC.h"
//MakeChairInstance2(tpcCorrection,St_TpcCurrentCorrectionC,Calibrations/tpc/TpcCurrentCorrection);
St_TpcCurrentCorrectionC *St_TpcCurrentCorrectionC::fgInstance = 0;
St_TpcCurrentCorrectionC *St_TpcCurrentCorrectionC::instance() {
  if (fgInstance) return fgInstance;
  St_tpcCorrection *table = (St_tpcCorrection *) StMaker::GetChain()->GetDataBase("Calibrations/tpc/TpcCurrentCorrectionX");
  if (! table)      table = (St_tpcCorrection *) StMaker::GetChain()->GetDataBase("Calibrations/tpc/TpcCurrentCorrection");
  assert(table);	  DEBUGTABLE(tpcCorrection);
  fgInstance = new St_TpcCurrentCorrectionC(table);
  return fgInstance;
}
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


St_MDFCorrectionC::St_MDFCorrectionC(St_MDFCorrection *table) : TChair(table), fFunc(0) {
  UInt_t N = table->GetNRows();
  fFunc = new TF1*[N];
  memset(fFunc, 0, N*sizeof(TF1*));
}


St_MDFCorrectionC::~St_MDFCorrectionC() {
  UInt_t N = Table()->GetNRows();
  for (UInt_t i = 0; i < N; i++) {SafeDelete(fFunc[i]);}
  delete [] fFunc;
}


Double_t St_MDFCorrectionC::MDFunc(Double_t *x, Double_t *p) {
  // Evaluate parameterization at point x. Optional argument coeff is
  // a vector of coefficients for the parameterisation, NCoefficients
  // elements long.
  assert(x);
  UInt_t k = p[0];
  assert(k >= 0 && k < fgMDFCorrectionC->getNumRows());
  Double_t returnValue = fgMDFCorrectionC->DMean(k);
  Double_t term        = 0;
  UChar_t    i, j;
  for (i = 0; i < fgMDFCorrectionC->NCoefficients(k); i++) {
    // Evaluate the ith term in the expansion
    term = fgMDFCorrectionC->Coefficients(k)[i];
    for (j = 0; j < fgMDFCorrectionC->NVariables(k); j++) {
      // Evaluate the factor (polynomial) in the j-th variable.
      Int_t    p  =  fgMDFCorrectionC->Powers(k)[i * fgMDFCorrectionC->NVariables(k) + j];
      Double_t y  =  1 + 2. / (fgMDFCorrectionC->XMax(k)[j] - fgMDFCorrectionC->XMin(k)[j])
	* (x[j] - fgMDFCorrectionC->XMax(k)[j]);
      term        *= fgMDFCorrectionC->EvalFactor(k,p,y);
    }
    // Add this term to the final result
    returnValue += term;
  }
  return returnValue;
}



Double_t St_MDFCorrectionC::Eval(Int_t k, Double_t x0, Double_t x1) const {
  Double_t x[2] = {x0, x1};
  return Eval(k,x);
}


Double_t St_MDFCorrectionC::Eval(Int_t k, Double_t *x) const {
  // Evaluate parameterization at point x. Optional argument coeff is
  // a vector of coefficients for the parameterisation, NCoefficients
  // elements long.
  assert(x);
  if (! fFunc[k]) {
    fgMDFCorrectionC = (St_MDFCorrectionC *) this;
    if (NVariables(k) <= 0) {
      return 0;
    } else if (NVariables(k) == 1) {
      fFunc[k] = new TF1(Form("%s_%i",Table()->GetName(),k),St_MDFCorrectionC::MDFunc,
			 XMin(k)[0],XMax(k)[0],1);
      fFunc[k]->SetParameter(0,k);
      fFunc[k]->Save(XMin(k)[0],XMax(k)[0],0,0,0,0);
    } else if (NVariables(k) == 2) {
      fFunc[k] = new TF2(Form("%s_%i",Table()->GetName(),k),St_MDFCorrectionC::MDFunc,
			 XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],1);
      fFunc[k]->SetParameter(0,k);
      ((TF2 *) fFunc[k])->Save(XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],0,0);
    } else if (NVariables(k) == 3) {
      fFunc[k] = new TF3(Form("%s_%i",Table()->GetName(),k),St_MDFCorrectionC::MDFunc,
			 XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],XMin(k)[2],XMax(k)[2],1);
      fFunc[k]->SetParameter(0,k);
      ((TF3 *) fFunc[k])->Save(XMin(k)[0],XMax(k)[0],XMin(k)[1],XMax(k)[1],XMin(k)[2],XMax(k)[2]);
    }
  }
  Double_t xx[3];
  for (Int_t v = 0; v < NVariables(k); v++) {
    xx[v] = TMath::Max(XMin(k)[v], TMath::Min(XMin(k)[v]+0.999*(XMax(k)[v]-XMin(k)[v]), x[v]));
  }
  Double_t returnValue = fFunc[k]->GetSave(xx);
  return returnValue;
}


Double_t St_MDFCorrectionC::EvalError(Int_t k, Double_t *x) const {
  // Evaluate parameterization error at point x. Optional argument coeff is
  // a vector of coefficients for the parameterisation, NCoefficients(k)
  // elements long.
  assert(x);
  Double_t returnValue = 0;
  Double_t term        = 0;
  UChar_t    i, j;
  for (i = 0; i < NCoefficients(k); i++) {
    // Evaluate the ith term in the expansion
    term = CoefficientsRMS(k)[i];
    for (j = 0; j < NVariables(k); j++) {
      // Evaluate the factor (polynomial) in the j-th variable.
      Int_t    p  =  Powers(k)[i * NVariables(k) + j];
      Double_t y  =  1 + 2. / (XMax(k)[j] - XMin(k)[j])
	* (x[j] - XMax(k)[j]);
      term        *= EvalFactor(p,y);
    }
    // Add this term to the final result
    returnValue += term*term;
  }
  returnValue = TMath::Sqrt(returnValue);
  return returnValue;
}


Double_t St_MDFCorrectionC::EvalFactor(Int_t k, Int_t p, Double_t x) const {
  // Evaluate function with power p at variable value x
  Int_t    i   = 0;
  Double_t p1  = 1;
  Double_t p2  = 0;
  Double_t p3  = 0;
  Double_t r   = 0;

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
#include "StDetectorDbMaker/St_tpcPedestalC.h"
MakeChairInstance(tpcPedestal,Calibrations/tpc/tpcPedestal);
#include "StDetectorDbMaker/St_tpcPadResponseC.h"
MakeChairInstance(tpcPadResponse,Calibrations/tpc/tpcPadResponse);
#include "StDetectorDbMaker/St_tpcSlowControlSimC.h"
MakeChairInstance(tpcSlowControlSim,Calibrations/tpc/tpcSlowControlSim);
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
tpcPadConfig_st *St_tpcPadConfigC::Struct(Int_t i)                        {return ((St_tpcPadConfig*) Table())->GetTable(i);}
UInt_t           St_tpcPadConfigC::getNumRows()                	          {return GetNRows();}
UChar_t          St_tpcPadConfigC::iTpc(Int_t sector)                     {UChar_t iTPC = Struct()->itpc[sector-1];  return iTPC;}
Int_t 	         St_tpcPadConfigC::padRows(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->padRows()               ;}
Int_t 	         St_tpcPadConfigC::innerPadRows(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->innerPadRows() 	   ;}
Int_t 	         St_tpcPadConfigC::innerPadRows48(Int_t sector) 	  {return St_tpcPadPlanesC::instance()->innerPadRows48()	   ;}
Int_t 	         St_tpcPadConfigC::innerPadRows52(Int_t sector) 	  {return St_tpcPadPlanesC::instance()->innerPadRows52()	   ;}
Int_t 	         St_tpcPadConfigC::outerPadRows(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->outerPadRows() 	   ;}
Int_t 	         St_tpcPadConfigC::superInnerPadRows(Int_t sector)        {return St_tpcPadPlanesC::instance()->superInnerPadRows()     ;}
Int_t 	         St_tpcPadConfigC::superOuterPadRows(Int_t sector)        {return St_tpcPadPlanesC::instance()->superOuterPadRows()     ;}
Double_t 	 St_tpcPadConfigC::innerSectorPadWidth(Int_t sector)      {return St_tpcPadPlanesC::instance()->innerSectorPadWidth()   ;}
Double_t 	 St_tpcPadConfigC::innerSectorPadLength(Int_t sector)     {return St_tpcPadPlanesC::instance()->innerSectorPadLength()  ;}
Double_t 	 St_tpcPadConfigC::innerSectorPadPitch(Int_t sector)      {return St_tpcPadPlanesC::instance()->innerSectorPadPitch()   ;}
Double_t 	 St_tpcPadConfigC::innerSectorRowPitch1(Int_t sector)     {return St_tpcPadPlanesC::instance()->innerSectorRowPitch1()  ;}
Double_t 	 St_tpcPadConfigC::innerSectorRowPitch2(Int_t sector)     {return St_tpcPadPlanesC::instance()->innerSectorRowPitch2()  ;}
Double_t 	 St_tpcPadConfigC::firstPadRow(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->firstPadRow() 	   ;}
Double_t 	 St_tpcPadConfigC::firstOuterSectorPadRow(Int_t sector)   {return St_tpcPadPlanesC::instance()->firstOuterSectorPadRow();}
Double_t 	 St_tpcPadConfigC::lastOuterSectorPadRow(Int_t sector)    {return St_tpcPadPlanesC::instance()->lastOuterSectorPadRow() ;}
Double_t 	 St_tpcPadConfigC::firstRowWidth(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->firstRowWidth()         ;}
Double_t 	 St_tpcPadConfigC::lastRowWidth(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->lastRowWidth()          ;}
Double_t 	 St_tpcPadConfigC::outerSectorPadWidth(Int_t sector)      {return St_tpcPadPlanesC::instance()->outerSectorPadWidth()   ;}
Double_t 	 St_tpcPadConfigC::outerSectorPadLength(Int_t sector)     {return St_tpcPadPlanesC::instance()->outerSectorPadLength()  ;}
Double_t 	 St_tpcPadConfigC::outerSectorPadPitch(Int_t sector)      {return St_tpcPadPlanesC::instance()->outerSectorPadPitch()   ;}
Double_t 	 St_tpcPadConfigC::outerSectorRowPitch(Int_t sector)      {return St_tpcPadPlanesC::instance()->outerSectorRowPitch()   ;}
Double_t 	 St_tpcPadConfigC::outerSectorLength(Int_t sector)        {return St_tpcPadPlanesC::instance()->outerSectorLength()     ;}
Double_t 	 St_tpcPadConfigC::ioSectorSeparation(Int_t sector)       {return St_tpcPadPlanesC::instance()->ioSectorSeparation()    ;}
Double_t 	 St_tpcPadConfigC::innerSectorEdge(Int_t sector) 	  {return St_tpcPadPlanesC::instance()->innerSectorEdge()       ;}
Double_t 	 St_tpcPadConfigC::outerSectorEdge(Int_t sector) 	  {return St_tpcPadPlanesC::instance()->outerSectorEdge() 	   ;}
Double_t 	 St_tpcPadConfigC::innerSectorPadPlaneZ(Int_t sector)     {return St_tpcPadPlanesC::instance()->innerSectorPadPlaneZ()  ;}
Double_t 	 St_tpcPadConfigC::outerSectorPadPlaneZ(Int_t sector)     {return St_tpcPadPlanesC::instance()->outerSectorPadPlaneZ()  ;}
Int_t* 	         St_tpcPadConfigC::innerPadsPerRow(Int_t sector) 	  {return St_tpcPadPlanesC::instance()->innerPadsPerRow()       ;}
Int_t* 	         St_tpcPadConfigC::outerPadsPerRow(Int_t sector) 	  {return St_tpcPadPlanesC::instance()->outerPadsPerRow()       ;}
Int_t            St_tpcPadConfigC::padsPerRow(Int_t sector, Int_t row)    {
  Int_t Ninner = innerPadRows(sector);
  return (row <= Ninner) ?
    innerPadsPerRow(sector)[row-1] :
    outerPadsPerRow(sector)[row-1-Ninner];
}
Double_t* 	 St_tpcPadConfigC::innerRowRadii(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->innerRowRadii()         ;}
Double_t* 	 St_tpcPadConfigC::outerRowRadii(Int_t sector) 	          {return St_tpcPadPlanesC::instance()->outerRowRadii() 	   ;}
// taken from StRItpcPadPlane
Int_t            St_tpcPadConfigC::numberOfRows(Int_t sector)             {return padRows(sector);}
Int_t            St_tpcPadConfigC::numberOfInnerRows(Int_t sector)        {return innerPadRows(sector);}
Int_t            St_tpcPadConfigC::numberOfInnerRows48(Int_t sector)      {return innerPadRows48(sector);}
Int_t            St_tpcPadConfigC::numberOfInnerRows52(Int_t sector)      {return innerPadRows52(sector);}
Int_t            St_tpcPadConfigC::numberOfOuterRows(Int_t sector)        {return outerPadRows(sector);}
Bool_t           St_tpcPadConfigC::isRowInRange(Int_t sector, Int_t row)  {return (row >= 1 && row<=numberOfRows(sector)) ? kTRUE: kFALSE;}
Double_t         St_tpcPadConfigC::radialDistanceAtRow(Int_t sector, Int_t row)       {
  if (! isRowInRange(sector,row)) return 0;
  Int_t Ninner = innerPadRows(sector);
  if ( row<=Ninner ) return innerRowRadii(sector)[row-1];
  else               return outerRowRadii(sector)[row-1-Ninner];
}
Int_t            St_tpcPadConfigC::numberOfPadsAtRow(Int_t sector, Int_t row)       {
  if (! isRowInRange(sector, row)) return 0;
  Int_t Ninner = innerPadRows(sector);
  if ( row<=Ninner ) return innerPadsPerRow(sector)[row-1];
  return outerPadsPerRow(sector)[row-1-Ninner];
}
Double_t         St_tpcPadConfigC::PadWidthAtRow(Int_t sector, Int_t row)       {
  if (! isRowInRange(sector,row)) return 0;
  Int_t Ninner = innerPadRows(sector);
  if ( row<=Ninner) return innerSectorPadWidth(sector);
  return outerSectorPadWidth(sector);
}
Double_t         St_tpcPadConfigC::PadLengthAtRow(Int_t sector, Int_t row)       {
  if (! isRowInRange(sector,row)) return 0;
  Int_t Ninner = innerPadRows(sector);
  if ( row<=Ninner) return innerSectorPadLength(sector);
  return outerSectorPadLength(sector);
}
Double_t         St_tpcPadConfigC::PadPitchAtRow(Int_t sector, Int_t row)       {
  if (! isRowInRange(sector,row)) return 0;
  Int_t Ninner = innerPadRows(sector);
  if ( row<=Ninner) return innerSectorPadPitch(sector);
  return outerSectorPadPitch(sector);
}
Double_t         St_tpcPadConfigC::RowPitchAtRow(Int_t sector, Int_t row)       {
  if (! isRowInRange(sector,row)) return 0;
  Int_t Ninner = innerPadRows(sector);
  if ( row<=numberOfInnerRows48(sector) ) return innerSectorRowPitch1(sector);
  else if (row>numberOfInnerRows48(sector)&&row<=Ninner) return innerSectorRowPitch2(sector);
  return outerSectorRowPitch(sector);
}
Int_t            St_tpcPadConfigC::indexForRowPad(Int_t sector, Int_t row, Int_t pad)       {
  if (pad >numberOfPadsAtRow(sector,row)) return -1;
  Int_t index = 0;
  Int_t Ninner = innerPadRows(sector);
  if (row>0 && row<=Ninner )             for (Int_t i=1;i<row;i++) index += numberOfPadsAtRow(sector,i);
  else
    if (row>Ninner&&row<=numberOfRows(sector)) for (Int_t i=Ninner+1;i<row;i++)  index += numberOfPadsAtRow(sector,i);
  index+=pad-1;
  return index;
}
#include "StDetectorDbMaker/St_TpcAvgPowerSupplyC.h"


void  St_tpcAnodeHVC::sockets(Int_t sector, Int_t padrow, Int_t &e1, Int_t &e2, Float_t &f2) {
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


Float_t St_tpcAnodeHVC::voltage(Int_t i) const {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
     LOG_ERROR << "St_tpcAnodeHVC::voltage(" << i << " is called but the valid St_TpcAvgPowerSupplyC::instance() exists\n";
  }
  return Struct(i)->voltage;
}


Float_t St_tpcAnodeHVC::voltagePadrow(Int_t sector, Int_t padrow) const {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->voltagePadrow(sector,padrow);
  }
  Int_t e1 = 0, e2 = 0;
  Float_t f2 = 0;
  St_tpcAnodeHVC::sockets(sector, padrow, e1, e2, f2);
  if (e1==0) return -99;
  Float_t v1=voltage(e1-1);
  if (f2 < 0.1) return v1;
  Float_t v2=voltage(e2-1);
  if (TMath::Abs(v2 - v1) > 40) return -99;
  if (TMath::Abs(v2 - v1) <  1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  Float_t B = (padrow <= St_tpcPadConfigC::instance()->innerPadRows(sector) ? 13.05e-3 : 10.26e-3);
  Float_t v_eff = TMath::Log((1.0-f2)*TMath::Exp(B*v1) + f2*TMath::Exp(B*v2)) / B;
  return v_eff;
}
MakeChairOptionalInstance(TpcAvgPowerSupply,Calibrations/tpc/TpcAvgPowerSupply);


Float_t St_TpcAvgPowerSupplyC::voltagePadrow(Int_t sector, Int_t padrow) const {
  Int_t e1 = 0, e2 = 0;
  Float_t f2 = 0;
  St_tpcAnodeHVC::sockets(sector, padrow, e1, e2, f2);
  if (e1==0) return -99;
  Int_t ch1 = St_TpcAvgCurrentC::ChannelFromSocket((e1-1)%19+1);
  Float_t v1=Voltage()[8*(sector-1)+ch1-1] ;
  if (f2==0) return v1;
  Int_t ch2 = St_TpcAvgCurrentC::ChannelFromSocket((e2-1)%19 + 1);
  if (ch1 == ch2) return v1;
  Float_t v2=Voltage()[8*(sector-1)+ch2-1] ;
  if (v2==v1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  Float_t B = (padrow <= St_tpcPadConfigC::instance()->innerPadRows(sector) ? 13.05e-3 : 10.26e-3);
  Float_t v_eff = TMath::Log((1.0-f2)*TMath::Exp(B*v1) + f2*TMath::Exp(B*v2)) / B;
  return v_eff;
}


Float_t St_TpcAvgPowerSupplyC::AcChargeL(Int_t sector, Int_t channel) {
  //  static const Double_t RA[2]        = { 154.484, 81.42}; // Outer/ Inner average Radii
  //  static const Double_t WireLenth[2] = {   3.6e5, 1.6e5};
  // L Inner = 190222, Outer = 347303
  static Float_t Length[8] = {
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


Float_t St_tpcAnodeHVavgC::voltage(Int_t i) const {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
     LOG_ERROR << "St_tpcAnodeHVavgC::voltage(" << i << " is called but the valid St_TpcAvgPowerSupplyC::instance() exists\n";
  }
  return Struct(i)->voltage;
}


Bool_t St_tpcAnodeHVavgC::tripped(Int_t sector, Int_t padrow) const {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->tripped(sector,padrow);
  }
  return (voltage() < -100);
}


Float_t St_tpcAnodeHVavgC::voltagePadrow(Int_t sector, Int_t padrow) const {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->voltagePadrow(sector,padrow);
  }
  Int_t e1 = 0, e2 = 0;
  Float_t f2 = 0;
  St_tpcAnodeHVC::sockets(sector, padrow, e1, e2, f2);
  if (e1==0) return -99;
  Float_t v1=voltage(e1-1);
  if (f2==0) return v1;
  Float_t v2=voltage(e2-1);
  if (v2==v1) return v1;
  // different voltages on influencing HVs
  // effective voltage is a sum of exponential gains
  Float_t B = (padrow <= St_tpcPadConfigC::instance()->innerPadRows(sector) ? 13.05e-3 : 10.26e-3);
  Float_t v_eff = TMath::Log((1.0-f2)*TMath::Exp(B*v1) + f2*TMath::Exp(B*v2)) / B;
  return v_eff;
}


#include "StDetectorDbMaker/St_tpcPadGainT0C.h"
MakeChairInstance(tpcPadGainT0,Calibrations/tpc/tpcPadGainT0);
#include "StDetectorDbMaker/St_itpcPadGainT0C.h"
MakeChairInstance(itpcPadGainT0,Calibrations/tpc/itpcPadGainT0);
#include "StDetectorDbMaker/St_tpcPadGainT0BC.h"
St_tpcPadGainT0BC *St_tpcPadGainT0BC::fgInstance = 0;
St_tpcPadGainT0BC *St_tpcPadGainT0BC::instance() {if (! fgInstance) fgInstance = new St_tpcPadGainT0BC(); return fgInstance;}


Float_t 	St_tpcPadGainT0BC::Gain(Int_t sector, Int_t row, Int_t pad) const {
  Float_t gain = 0;
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


Float_t 	  St_tpcPadGainT0BC::T0(Int_t sector, Int_t row, Int_t pad) const {
  Float_t T0 = 0;
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


Bool_t    St_tpcPadGainT0BC::livePadrow(Int_t sector, Int_t row) const {
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


Int_t St_TpcAvgCurrentC::ChannelFromRow(Int_t sector, Int_t row) {
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


Int_t St_TpcAvgCurrentC::ChannelFromSocket(Int_t socket) {
  Int_t channel = -1;
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


Float_t St_TpcAvgCurrentC::AcChargeL(Int_t sector, Int_t channel) {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->AcChargeL(sector,channel);
  }
  //  static const Double_t RA[2]        = { 154.484, 81.42}; // Outer/ Inner average Radii
  //  static const Double_t WireLenth[2] = {   3.6e5, 1.6e5};
  // L Inner = 190222, Outer = 347303
  static Float_t Length[8] = {
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


Float_t St_TpcAvgCurrentC::AvCurrent(Int_t sector, Int_t channel) {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->AvCurrent(sector,channel);
  }
  return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
    Struct()->AvCurrent[8*(sector-1)+channel-1] :     0;
}


Float_t St_TpcAvgCurrentC::AcCharge(Int_t sector, Int_t channel) {
  if (! St_TpcAvgPowerSupplyC::instance()->Table()->IsMarked()) {
    return St_TpcAvgPowerSupplyC::instance()->AcCharge(sector,channel);
  }
  return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
    Struct()->AcCharge[8*(sector-1)+channel-1] :     0;
}
#include "StDetectorDbMaker/St_tpcRDOMapC.h"
MakeChairInstance(tpcRDOMap,Calibrations/tpc/tpcRDOMap);


Int_t St_tpcRDOMapC::rdo(Int_t padrow, Int_t pad) const {
  Int_t rdo = 0;
  Int_t N = nrows(0);
  for (Int_t i = 0; i < N; i++) {
    if (padrow != row(i)) continue;
    if (pad < padMin(i) || pad > padMax(i)) continue;
    rdo = rdoI(i);
    
    break;
  }
  return rdo;
}
#include "StDetectorDbMaker/St_tpcRDOT0offsetC.h"
MakeChairInstance(tpcRDOT0offset,Calibrations/tpc/tpcRDOT0offset);
Float_t St_tpcRDOT0offsetC::T0(Int_t sector, Int_t padrow, Int_t pad) const {
  Float_t t0 = 0;
  if (! IsShfited(sector)) return t0;
  if (St_tpcPadConfigC::instance()->iTPC(sector) && padrow <= 40)  return t0; // no shift in iTPC
  Int_t rdo = St_tpcRDOMapC::instance()->rdo(padrow,pad);
  if (!rdo) return t0;
  t0 = Struct()->t0[sector-1][rdo-1];
  return t0;
}

#include "StDetectorDbMaker/St_trigDetSumsC.h"
MakeChairInstance(trigDetSums, Calibrations/rich/trigDetSums);
St_trigDetSums* St_trigDetSumsC::fgTableCopy = 0;
//___________________tpc_____________________________________________________________
#include "StDetectorDbMaker/St_tss_tssparC.h"
MakeChairInstance(tss_tsspar,tpc/tsspars/tsspar);


Float_t St_tss_tssparC::gain(Int_t sector, Int_t row) {
  Int_t l = 0;
  Double_t V_nominal = 1390;
  Float_t V = 0;
  Float_t gain = 0;
  if (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) {l = 1; V_nominal = 1170;}
  St_tpcGainCorrectionC *gC = St_tpcGainCorrectionC::instance();
  Int_t NRows = gC->GetNRows();
  if (l >= NRows) return gain;
  V = St_tpcAnodeHVavgC::instance()->voltagePadrow(sector,row);
  if (V > 0) {
    Double_t v = V - V_nominal;
#if 0
    // Hack for Run XVI
    if ( gC->min(l) > -1 && v > -100 && v < 0) {
      if (l == 0) v =   0;
      else        v = -70;
    } else
#endif
    if (v < gC->min(l) || v > gC->max(l)) return gain;
    if (gC->min(l) < -450) {
      // if range was expanded below 150 V then use only the linear approximation
      gain  = TMath::Exp(gC->CalcCorrection(l,v, 0., 2));
    } else {
      gain  = TMath::Exp(gC->CalcCorrection(l,v));
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


starClockOnl_st *St_starClockOnlC::Struct(Int_t i) {
  starClockOnl_st *s = ((St_starClockOnl* ) instance()->Table())->GetTable();
  Int_t N =  getNumRows(); // with i < 0 look for positive frequency
  if (i >= 0 && i < N) return s + i;
  for (Int_t j = 0; j < N; j++, s++) if (s->frequency > 0) break;
  assert(s->frequency > 0 && s->frequency < 1e7);
  return s;
}
#include "StDetectorDbMaker/St_starMagOnlC.h"
MakeChairInstance(starMagOnl,RunLog/onl/starMagOnl);
#include "St_beamInfoC.h"
MakeChairInstance(beamInfo,RunLog/onl/beamInfo);


Bool_t        St_beamInfoC::IsFixedTarget() {
  Bool_t isFixTag = kFALSE;
  Float_t MaxIntensity = TMath::Max(blueIntensity(), yellowIntensity());
  Float_t MinIntensity = TMath::Min(blueIntensity(), yellowIntensity());
  if (MaxIntensity > 1.0 && MaxIntensity > 10*MinIntensity) isFixTag = kTRUE;
  return isFixTag;
}


#include "StDetectorDbMaker/St_tpcRDOMasksC.h"
MakeChairInstance(tpcRDOMasks,RunLog/onl/tpcRDOMasks);


UInt_t       St_tpcRDOMasksC::getSectorMask(UInt_t sector) {
  UInt_t MASK = 0x0000; // default is to mask it out
  //UInt_t MASK = 0xFFFF; // change to  ON by default ** THIS WAS A HACK
  if(sector < 1 || sector > 24 || getNumRows() == 0){
    LOG_WARN << "St_tpcRDOMasksC:: getSectorMask : return default mask for "
     << "sector= " << sector << " getNumRows()=" << getNumRows() << '\n';
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
#include "StDetectorDbMaker/St_tpcSectorPositionC.h"
St_tpcSectorPositionC *St_tpcSectorPositionC::fgInstance = 0;
St_tpcSectorPosition  *St_tpcSectorPositionC::fgTables[24] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
St_tpcSectorPositionC *St_tpcSectorPositionC::instance() {
  if (fgInstance) return fgInstance;
  TDataSet *tpc = StMaker::GetChain()->GetDataBase("Geometry/tpc");
  assert(tpc);
  for (Int_t sec = 1; sec <= 24; sec++) {
    TString path = Form("Sector_%02i/tpcSectorPosition",sec);
    fgTables[sec-1] = (St_tpcSectorPosition  *) tpc->Find(path.Data());
    if (! fgTables[sec-1]) {
      LOG_WARN << "St_tpcSectorPositionC::instance " << tpc->GetTitle() << "/" << path.Data()
       << "\twas not found" << endm;
      assert(fgTables[sec-1]);
    }
    {
      TDatime t[2];
      St_db_Maker::GetValidity(fgTables[sec-1],t);
      Int_t Nrows = fgTables[sec-1]->GetNRows();
      LOG_WARN << "St_tpcSectorPositionC::instance found table " << fgTables[sec-1]->GetName()
	       << " with NRows = " << Nrows << " in db" << endm;
      LOG_WARN << "Validity:" << t[0].GetDate() << "/" << t[0].GetTime()
	       << " -----   " << t[1].GetDate() << "/" << t[1].GetTime() << endm;
      fgTables[sec-1]->Print(0,1);
    }
  }
  fgInstance = new St_tpcSectorPositionC();
  return fgInstance;
}
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

St_SurveyC::St_SurveyC(St_Survey *table) : TChair(table), fRotations(0)  {
  UInt_t N = getNumRows();
  fRotations = new TGeoHMatrix*[N];
  for (UInt_t i = 0; i < N; i++) {
    fRotations[i] = new TGeoHMatrix;
    TGeoHMatrix &rot = *fRotations[i];
    if (N == 1) rot.SetName(Table()->GetName());
    else        rot.SetName(Form("%s_%i",Table()->GetName(),i+1));
    rot.SetRotation(Rotation(i));
    rot.SetTranslation(Translation(i));
    Normalize(rot);
    assert(TMath::Abs(rot.Determinant())-1 < 1.e-3);
  }
}


St_SurveyC::~St_SurveyC() {
  if (fRotations) {
    for (UInt_t i = 0; i < getNumRows(); i++) {
      SafeDelete(fRotations[0]);
    }
    SafeDelete(fRotations);
  }
}


Double_t St_SurveyC::IsOrtogonal(const Double_t *r) {
// Perform orthogonality test for rotation.
  Double_t cmax = 0;
  Double_t cij;
  for (Int_t i=0; i<2; i++) {
    for (Int_t j=i+1; j<3; j++) {
      // check columns
      cij = TMath::Abs(r[i]*r[j]+r[i+3]*r[j+3]+r[i+6]*r[j+6]);
      if (cij>1E-4) cmax = cij;
      // check rows
      cij = TMath::Abs(r[3*i]*r[3*j]+r[3*i+1]*r[3*j+1]+r[3*i+2]*r[3*j+2]);
      if (cij>cmax) cmax = cij;
    }
  }
  return cmax;
}


void St_SurveyC::Normalize(TGeoHMatrix &R) {
#if 0
  Double_t det = R.Determinant();
  Double_t ort = IsOrtogonal(R.GetRotationMatrix());
  static Double_t eps = 1e-7;
  if ( TMath::Abs(TMath::Abs(det) - 1) < eps && ort < eps) return;
  LOG_INFO << "St_SurveyC::Normalize matrix " << R.GetName()
   << Form(" has determinant-1 = %10.7f\tortoganality %10.7f",TMath::Abs(det)-1,ort) << endm;
  cout << "Old\t"; R.Print();
  const Double_t *r = R.GetRotationMatrix();
  SMatrix<double,3,3> A(r,9); //   cout << "A: " << endl << A << endl;
  SMatrix<double,3,3> B = A;
  A.Det(det); //  cout << "Determinant - 1: " << det-1 << endl;  cout << "A again: " << endl << A << endl;
  A = B;
  A.Invert();//   cout << "A^-1: " << endl << A << endl;
  // check if this is really the inverse:  cout << "A^-1 * B: " << endl << A * B << endl;
  // the Babylonian method for extracting the square root of a matrix :  Q_{n+1} = 2 * M * ((Q_{n}^{-1} * M) + (M^{T} *Q_{n}))^{-1}
  SMatrix<double,3,3> Qn1;
  SMatrix<double,3,3> Qn2;
  SMatrix<double,3,3> M = B;
  SMatrix<double,3,3> Qn = M;
  Int_t ifail = 0;
  Int_t N = 0;
  Qn.Det(det); if (_debug) {LOG_INFO << "N " << N << "\tQn Determinant - 1: " << Form("%15.5g",det-1) << endm;}
  Qn = M;
  while (TMath::Abs(TMath::Abs(det) - 1) > eps) {
    SMatrix<double,3> QnInv = Qn.Inverse(ifail);
    if (ifail) {
      LOG_ERROR << "St_SurveyC::Normalize:: Qn inversion failed" << endm;
      break;
    }
    SMatrix<double,3,3> C1 = QnInv * M;
    SMatrix<double,3,3> C2 = Transpose(M) * Qn;
    SMatrix<double,3,3> C  = C1 + C2;
    SMatrix<double,3,3> CInv = C.Inverse(ifail);
    if (ifail) {
      LOG_ERROR << "St_SurveyC::Normalize:: C inversion failed" << endm;
      break;
    }
    Qn1 = 2 * M * CInv;
    Qn2 = Qn1;
    N++;
    Qn2.Det(det);  if (_debug) {LOG_INFO << "N " << N << "\tQn2 Determinant - 1: " << Form("%15.5g",det-1) << endm;}
    if (N > 13) break;
    Qn = Qn1;
    if (_debug) {LOG_INFO << "Qn:" << endl << Qn << endm;}
  }
  R.SetRotation(Qn.Array()); cout << "New\t"; R.Print();
#endif
  if (_debug) {
    LOG_INFO << "Matrix:\n"; R.Print("");
    LOG_INFO << "Determinant-1 = " << R.Determinant()-1 << '\n';
    const Double_t *rr = R.GetRotationMatrix();
    LOG_INFO << "Ortogonality " << IsOrtogonal(rr) << '\n';
  }
  return;
}


const TGeoHMatrix &St_SurveyC::GetMatrix(Int_t i) {
  assert(fRotations || fRotations[i]);
  assert(TMath::Abs(fRotations[i]->Determinant())-1 < 1.e-3);
  return *fRotations[i];
}


const TGeoHMatrix &St_SurveyC::GetMatrix4Id(Int_t id) {
  for (UInt_t i = 0; i < getNumRows(); i++) {
    if (Id(i) == id) {
      return GetMatrix(i);
    }
  }
  LOG_INFO  << "St_SurveyC::GetMatrix4Id(" << id << ") entry has not been found\n";
  const TTable *table = Table();
  Int_t Nrows = table->GetNRows();
  table->Print(0,Nrows);
  assert(0);
  return GetMatrix(0);
}


const TGeoHMatrix &St_SurveyC::GetMatrixR(Int_t i) {
  static TGeoHMatrix rot;
  Double_t rotations[9] = {
    r00(i), r01(i),      0,
    r10(i), r11(i),      0,
    0     ,      0, r22(i)};
  rot.SetName(Form("%s_%i",Table()->GetName(),i));
  rot.SetRotation(rotations);
  rot.SetTranslation(Translation(i));
  return *&rot;
}


void St_SurveyC::GetAngles(Double_t &phi, Double_t &the, Double_t &psi, Int_t i) {
  phi = the = psi = 0;  // Korn 14.10-5
  Double_t cosDelta = (r00(i) + r11(i) + r22(i) - 1)/2; // (Tr(R) - 1)/2
  Double_t Delta = TMath::ACos(cosDelta);
  if (Delta < 0) Delta += 2*TMath::Pi();
  Double_t sinDelta2 = TMath::Sin(Delta/2);
  if (TMath::Abs(sinDelta2) < 1.e-7) return;
  Double_t c[3] = {
    (r21(i) - r12(i))/(2*sinDelta2), // a32-a23
    (r02(i) - r20(i))/(2*sinDelta2), // a13-a31
    (r10(i) - r01(i))/(2*sinDelta2)  // a21-a12
  };
  Double_t u = TMath::ATan2(c[0],c[1]);
  Double_t v = TMath::ATan(c[2]*TMath::Tan(Delta/2));
  phi = (v - u)/2 - TMath::Pi()/2;
  psi = (v + u)/2 - TMath::Pi()/2;
  the = 2*TMath::ATan2(c[0]*TMath::Sin(v),c[2]*TMath::Sin(u));
  Double_t raddeg = 180./TMath::Pi();
  phi   *= raddeg;
  the   *= raddeg;
  psi   *= raddeg;
}
