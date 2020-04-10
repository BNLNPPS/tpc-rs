#ifndef St_richvoltagesC_h
#define St_richvoltagesC_h

#include "tpcrs/config_structs.h"
#include "richvoltages.h"

struct St_richvoltagesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_richvoltagesC, richvoltages_st>
{
  UInt_t 	runNumber(Int_t i = 0) 	        {return Struct(i)->runNumber;}
  UInt_t 	startStatusTime(Int_t i = 0) 	{return Struct(i)->startStatusTime;}
  UInt_t 	endStatusTime(Int_t i = 0) 	{return Struct(i)->endStatusTime;}
  UInt_t 	status(Int_t i = 0) 	        {return Struct(i)->status;}
};
#endif
