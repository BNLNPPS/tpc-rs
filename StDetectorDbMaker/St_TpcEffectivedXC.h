#ifndef St_TpcEffectivedXC_h
#define St_TpcEffectivedXC_h

#include "tpcrs/config_structs.h"
#include "TpcEffectivedX.h"

struct St_TpcEffectivedXC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcEffectivedXC, TpcEffectivedX_st> {
  float 	scaleInner(int i = 0) 	const {return Struct(i)->scaleInner;}
  float 	scaleOuter(int i = 0) 	const {return Struct(i)->scaleOuter;}
};
#endif
