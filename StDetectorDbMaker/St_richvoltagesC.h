#ifndef St_richvoltagesC_h
#define St_richvoltagesC_h

#include "tpcrs/config_structs.h"
#include "richvoltages.h"

struct St_richvoltagesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_richvoltagesC, richvoltages_st>
{
  unsigned int 	runNumber(int i = 0) 	        {return Struct(i)->runNumber;}
  unsigned int 	startStatusTime(int i = 0) 	{return Struct(i)->startStatusTime;}
  unsigned int 	endStatusTime(int i = 0) 	{return Struct(i)->endStatusTime;}
  unsigned int 	status(int i = 0) 	        {return Struct(i)->status;}
};
#endif
