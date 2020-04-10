#ifndef St_MDFCorrectionC_h
#define St_MDFCorrectionC_h
#include <string>
#include "tpcrs/config_structs.h"
#include "MDFCorrection.h"
#include "TF3.h"
#include "TF2.h"
#include "TF1.h"
struct St_MDFCorrectionC : tpcrs::IConfigStruct {
  virtual MDFCorrection_st* Struct(int i = 0) const = 0;
  enum EMDFPolyType {
    kMonomials,
    kChebyshev,
    kLegendre
  };
  St_MDFCorrectionC();

  void Initialize()
  {
    UInt_t N = GetNRows();
    fFunc = new TF1*[N];
    memset(fFunc, 0, N*sizeof(TF1*));
  }

  UChar_t 	idx(Int_t k = 0)        	const {return Struct(k)->idx;}
  UChar_t 	nrows(Int_t k = 0) 	        const {return Struct(k)->nrows;}
  UChar_t 	PolyType(Int_t k = 0) 	        const {return Struct(k)->PolyType;}
  UChar_t 	NVariables(Int_t k = 0) 	const {return Struct(k)->NVariables;}
  UChar_t 	NCoefficients(Int_t k = 0) 	const {return Struct(k)->NCoefficients;}
  UChar_t* 	Powers(Int_t k = 0) 	        const {return Struct(k)->Power;}
  Double_t 	DMean(Int_t k = 0)           	const {return Struct(k)->DMean;}
  Double_t* 	XMin(Int_t k = 0)         	const {return Struct(k)->XMin;}
  Double_t* 	XMax(Int_t k = 0)       	const {return Struct(k)->XMax;}
  Double_t* 	Coefficients(Int_t k = 0) 	const {return Struct(k)->Coefficients;}
  Double_t* 	CoefficientsRMS(Int_t k = 0) 	const {return Struct(k)->CoefficientsRMS;}
  Double_t      Eval(Int_t k = 0, Double_t *x = 0) const;
  Double_t      Eval(Int_t k, Double_t x0, Double_t x1) const;
  Double_t      EvalError(Int_t k = 0, Double_t *x = 0) const;
  static Double_t MDFunc(Double_t *x = 0, Double_t *p = 0);
  static St_MDFCorrectionC *fgMDFCorrectionC;
 protected:
  virtual ~St_MDFCorrectionC();
 private:
  Double_t EvalFactor(Int_t k = 0, Int_t p = 0, Double_t x = 0) const;
  TF1         **fFunc;
};
#endif
