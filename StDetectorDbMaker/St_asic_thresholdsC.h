#ifndef St_asic_thresholdsC_h
#define St_asic_thresholdsC_h

#include "tpcrs/config_structs.h"
#include "asic_thresholds.h"

struct St_asic_thresholdsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_asic_thresholdsC, asic_thresholds_st>
{
  int 	thresh_lo(int i = 0) 	{return Struct(i)->thresh_lo;}
  int 	thresh_hi(int i = 0) 	{return Struct(i)->thresh_hi;}
  int 	n_seq_lo(int i = 0) 	{return Struct(i)->n_seq_lo;}
  int 	n_seq_hi(int i = 0) 	{return Struct(i)->n_seq_hi;}
};
#endif
