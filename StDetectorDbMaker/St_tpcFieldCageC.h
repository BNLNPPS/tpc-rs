#ifndef St_tpcFieldCageC_h
#define St_tpcFieldCageC_h

#include "tpcrs/config_structs.h"
#include "tpcFieldCage.h"

struct St_tpcFieldCageC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcFieldCageC, tpcFieldCage_st>
{
  float 	innerFieldCageShift(int i = 0) {return Struct(i)->innerFieldCageShift;}
  float 	InnerFieldCageShift(int i = 0) {return innerFieldCageShift(i);}
  float 	eastClockError(int i = 0) 	{return Struct(i)->eastClockError;}
  float 	EastClockError(int i = 0) 	{return eastClockError(i);}
  float 	westClockError(int i = 0) 	{return Struct(i)->westClockError;}
  float 	WestClockError(int i = 0) 	{return westClockError(i);}
};
#endif
