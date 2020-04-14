#ifndef St_tpcFieldCageShortC_h
#define St_tpcFieldCageShortC_h

#include "tpcrs/config_structs.h"
#include "tpcFieldCageShort.h"

struct St_tpcFieldCageShortC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcFieldCageShortC, tpcFieldCageShort_st> {
  float 	side(int i = 0) 	        {return Struct(i)->side;}
  float 	cage(int i = 0) 	        {return Struct(i)->cage;}
  float 	ring(int i = 0) 	        {return Struct(i)->ring;}
  float 	resistor(int i = 0) 	        {return Struct(i)->resistor;}
  float 	MissingResistance(int i = 0) 	{return Struct(i)->MissingResistance;}
};
#endif
