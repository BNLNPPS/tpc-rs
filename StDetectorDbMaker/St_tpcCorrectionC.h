#ifndef St_tpcCorrectionC_h
#define St_tpcCorrectionC_h

#include "tpcrs/config_structs.h"
#include "tpcCorrection.h"

struct St_tpcCorrectionC : tpcrs::IConfigStruct
{
  virtual tpcCorrection_st* Struct(int i = 0) const = 0;
  Int_t 	type(Int_t i = 0) 	const {return Struct(i)->type;}
  Int_t 	idx(Int_t i = 0) 	const {return Struct(i)->idx;}
  Int_t 	nrows(Int_t i = 0) 	const {return Struct(i)->nrows;}
  Int_t 	npar(Int_t i = 0) 	const {return Struct(i)->npar;}
  Double_t 	OffSet(Int_t i = 0) 	const {return Struct(i)->OffSet;}
  Double_t 	min(Int_t i = 0) 	const {return Struct(i)->min;}
  Double_t 	max(Int_t i = 0) 	const {return Struct(i)->max;}
  Double_t* 	a(Int_t i = 0) 	        const {return Struct(i)->a;}
  Double_t CalcCorrection(Int_t i, Double_t x, Double_t z = 0, Int_t NparMax = -1);
  Double_t SumSeries(tpcCorrection_st* cor, Double_t x, Double_t z = 0, Int_t NparMax = -1);
};
#endif
