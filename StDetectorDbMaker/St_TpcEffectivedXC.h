#ifndef St_TpcEffectivedXC_h
#define St_TpcEffectivedXC_h

#include "tpcrs/config_structs.h"
#include "TpcEffectivedX.h"

struct St_TpcEffectivedXC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcEffectivedXC, TpcEffectivedX_st> {
  Float_t 	scaleInner(Int_t i = 0) 	const {return Struct(i)->scaleInner;}
  Float_t 	scaleOuter(Int_t i = 0) 	const {return Struct(i)->scaleOuter;}
};
#endif
