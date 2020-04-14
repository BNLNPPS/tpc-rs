#ifndef St_tpcCorrectionC_h
#define St_tpcCorrectionC_h

#include "tpcrs/config_structs.h"
#include "tpcCorrection.h"

struct St_tpcCorrectionC : tpcrs::IConfigStruct
{
  virtual tpcCorrection_st* Struct(int i = 0) const = 0;
  int 	type(int i = 0) 	const {return Struct(i)->type;}
  int 	idx(int i = 0) 	const {return Struct(i)->idx;}
  int 	nrows(int i = 0) 	const {return Struct(i)->nrows;}
  int 	npar(int i = 0) 	const {return Struct(i)->npar;}
  double 	OffSet(int i = 0) 	const {return Struct(i)->OffSet;}
  double 	min(int i = 0) 	const {return Struct(i)->min;}
  double 	max(int i = 0) 	const {return Struct(i)->max;}
  double* 	a(int i = 0) 	        const {return Struct(i)->a;}
  double CalcCorrection(int i, double x, double z = 0, int NparMax = -1);
  double SumSeries(tpcCorrection_st* cor, double x, double z = 0, int NparMax = -1);
};
#endif
