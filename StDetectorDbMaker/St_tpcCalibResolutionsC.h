#ifndef St_tpcCalibResolutionsC_h
#define St_tpcCalibResolutionsC_h

#include "tpcrs/config_structs.h"
#include "tpcCalibResolutions.h"

struct St_tpcCalibResolutionsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcCalibResolutionsC, tpcCalibResolutions_st> {
  Float_t 	SpaceCharge(Int_t i = 0) 	{return Struct(i)->SpaceCharge;}
  Float_t 	GridLeak(Int_t i = 0)	 	{return Struct(i)->GridLeak;}
};
#endif
