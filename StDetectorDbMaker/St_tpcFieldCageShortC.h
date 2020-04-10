#ifndef St_tpcFieldCageShortC_h
#define St_tpcFieldCageShortC_h

#include "tpcrs/config_structs.h"
#include "tpcFieldCageShort.h"

struct St_tpcFieldCageShortC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcFieldCageShortC, tpcFieldCageShort_st> {
  Float_t 	side(Int_t i = 0) 	        {return Struct(i)->side;}
  Float_t 	cage(Int_t i = 0) 	        {return Struct(i)->cage;}
  Float_t 	ring(Int_t i = 0) 	        {return Struct(i)->ring;}
  Float_t 	resistor(Int_t i = 0) 	        {return Struct(i)->resistor;}
  Float_t 	MissingResistance(Int_t i = 0) 	{return Struct(i)->MissingResistance;}
};
#endif
