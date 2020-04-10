#ifndef St_tpcFieldCageC_h
#define St_tpcFieldCageC_h

#include "tpcrs/config_structs.h"
#include "tpcFieldCage.h"

struct St_tpcFieldCageC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcFieldCageC, tpcFieldCage_st>
{
  Float_t 	innerFieldCageShift(Int_t i = 0) {return Struct(i)->innerFieldCageShift;}
  Float_t 	InnerFieldCageShift(Int_t i = 0) {return innerFieldCageShift(i);}
  Float_t 	eastClockError(Int_t i = 0) 	{return Struct(i)->eastClockError;}
  Float_t 	EastClockError(Int_t i = 0) 	{return eastClockError(i);}
  Float_t 	westClockError(Int_t i = 0) 	{return Struct(i)->westClockError;}
  Float_t 	WestClockError(Int_t i = 0) 	{return westClockError(i);}
};
#endif
