#ifndef St_tpcDriftVelocityC_h
#define St_tpcDriftVelocityC_h

#include "tpcrs/config_structs.h"
#include "tpcDriftVelocity.h"

struct St_tpcDriftVelocityC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcDriftVelocityC, tpcDriftVelocity_st> {
  Float_t 	laserDriftVelocityEast(Int_t i = 0) 	{return Struct(i)->laserDriftVelocityEast;}
  Float_t 	laserDriftVelocityWest(Int_t i = 0) 	{return Struct(i)->laserDriftVelocityWest;}
  Float_t 	cathodeDriftVelocityEast(Int_t i = 0) 	{return Struct(i)->cathodeDriftVelocityEast;}
  Float_t 	cathodeDriftVelocityWest(Int_t i = 0) 	{return Struct(i)->cathodeDriftVelocityWest;}
};
#endif
