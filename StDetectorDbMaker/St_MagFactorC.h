#ifndef St_MagFactorC_h
#define St_MagFactorC_h

#include "tpcrs/config_structs.h"
#include "MagFactor.h"

struct St_MagFactorC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_MagFactorC, MagFactor_st> {
  float 	ScaleFactor(int i = 0) {return Struct(i)->ScaleFactor;}
};
#endif
