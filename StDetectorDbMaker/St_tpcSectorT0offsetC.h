#ifndef St_tpcSectorT0offsetC_h
#define St_tpcSectorT0offsetC_h

#include "tpcrs/config_structs.h"
#include "tpcSectorT0offset.h"

struct St_tpcSectorT0offsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcSectorT0offsetC, tpcSectorT0offset_st> {
  Float_t* 	t0(Int_t i = 0) 	        {return Struct(i)->t0;}
  Float_t       t0offset(Int_t sector=1)        {return t0()[sector-1];}
};
#endif
