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
    unsigned int N = GetNRows();
    fFunc = new TF1*[N];
    memset(fFunc, 0, N*sizeof(TF1*));
  }

  unsigned char 	idx(int k = 0)        	const {return Struct(k)->idx;}
  unsigned char 	nrows(int k = 0) 	        const {return Struct(k)->nrows;}
  unsigned char 	PolyType(int k = 0) 	        const {return Struct(k)->PolyType;}
  unsigned char 	NVariables(int k = 0) 	const {return Struct(k)->NVariables;}
  unsigned char 	NCoefficients(int k = 0) 	const {return Struct(k)->NCoefficients;}
  unsigned char* 	Powers(int k = 0) 	        const {return Struct(k)->Power;}
  double 	DMean(int k = 0)           	const {return Struct(k)->DMean;}
  double* 	XMin(int k = 0)         	const {return Struct(k)->XMin;}
  double* 	XMax(int k = 0)       	const {return Struct(k)->XMax;}
  double* 	Coefficients(int k = 0) 	const {return Struct(k)->Coefficients;}
  double* 	CoefficientsRMS(int k = 0) 	const {return Struct(k)->CoefficientsRMS;}
  double      Eval(int k = 0, double *x = 0) const;
  double      Eval(int k, double x0, double x1) const;
  double      EvalError(int k = 0, double *x = 0) const;
  static double MDFunc(double *x = 0, double *p = 0);
  static St_MDFCorrectionC *fgMDFCorrectionC;
 protected:
  virtual ~St_MDFCorrectionC();
 private:
  double EvalFactor(int k = 0, int p = 0, double x = 0) const;
  TF1         **fFunc;
};
#endif
