#ifndef St_TpcSecRowCorC_h
#define St_TpcSecRowCorC_h

#include "tpcrs/config_structs.h"
#include "TpcSecRowCor.h"

struct St_TpcSecRowCorC : tpcrs::IConfigStruct
{
  virtual TpcSecRowCor_st* Struct(int i = 0) const = 0;

  Float_t* 	GainScale(Int_t i = 0) 	        {return Struct(i)->GainScale;}
  Float_t* 	GainRms(Int_t i = 0) 	        {return Struct(i)->GainRms;}
};
#endif
