#ifndef St_tpcSectorT0offsetC_h
#define St_tpcSectorT0offsetC_h

#include "tpcrs/config_structs.h"
#include "tpcSectorT0offset.h"

struct St_tpcSectorT0offsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcSectorT0offsetC, tpcSectorT0offset_st> {
  float* 	t0(int i = 0) 	        {return Struct(i)->t0;}
  float       t0offset(int sector=1)        {return t0()[sector-1];}
};
#endif
