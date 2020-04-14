#ifndef St_tpcDriftVelocityC_h
#define St_tpcDriftVelocityC_h

#include "tpcrs/config_structs.h"
#include "tpcDriftVelocity.h"

struct St_tpcDriftVelocityC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcDriftVelocityC, tpcDriftVelocity_st> {
  float 	laserDriftVelocityEast(int i = 0) 	{return Struct(i)->laserDriftVelocityEast;}
  float 	laserDriftVelocityWest(int i = 0) 	{return Struct(i)->laserDriftVelocityWest;}
  float 	cathodeDriftVelocityEast(int i = 0) 	{return Struct(i)->cathodeDriftVelocityEast;}
  float 	cathodeDriftVelocityWest(int i = 0) 	{return Struct(i)->cathodeDriftVelocityWest;}
};
#endif
