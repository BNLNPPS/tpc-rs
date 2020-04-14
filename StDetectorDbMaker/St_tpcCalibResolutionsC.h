#ifndef St_tpcCalibResolutionsC_h
#define St_tpcCalibResolutionsC_h

#include "tpcrs/config_structs.h"
#include "tpcCalibResolutions.h"

struct St_tpcCalibResolutionsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcCalibResolutionsC, tpcCalibResolutions_st> {
  float 	SpaceCharge(int i = 0) 	{return Struct(i)->SpaceCharge;}
  float 	GridLeak(int i = 0)	 	{return Struct(i)->GridLeak;}
};
#endif
