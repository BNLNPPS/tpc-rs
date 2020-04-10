#ifndef St_asic_thresholdsC_h
#define St_asic_thresholdsC_h

#include "tpcrs/config_structs.h"
#include "asic_thresholds.h"

struct St_asic_thresholdsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_asic_thresholdsC, asic_thresholds_st>
{
  Int_t 	thresh_lo(Int_t i = 0) 	{return Struct(i)->thresh_lo;}
  Int_t 	thresh_hi(Int_t i = 0) 	{return Struct(i)->thresh_hi;}
  Int_t 	n_seq_lo(Int_t i = 0) 	{return Struct(i)->n_seq_lo;}
  Int_t 	n_seq_hi(Int_t i = 0) 	{return Struct(i)->n_seq_hi;}
};
#endif
