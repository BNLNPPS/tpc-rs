#ifndef St_tpcPadrowT0C_h
#define St_tpcPadrowT0C_h
#include "tpcrs/config_structs.h"
#include "tpcPadrowT0.h"

struct St_tpcPadrowT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadrowT0C, tpcPadrowT0_st> {
  Float_t T0(Int_t sector, Int_t row) {return Struct(sector-1)->T0[row-1];}
};

#endif
