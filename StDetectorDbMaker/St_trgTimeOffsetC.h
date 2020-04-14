#ifndef St_trgTimeOffsetC_h
#define St_trgTimeOffsetC_h

#include "tpcrs/config_structs.h"
#include "trgTimeOffset.h"

struct St_trgTimeOffsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_trgTimeOffsetC, trgTimeOffset_st>
{
  float 	offset(int i = 0)     	   {return Struct(i)->offset;}
  float 	laserOffset(int i = 0) 	   {return Struct(i)->laserOffset;}
  float 	laserOffsetW(int i = 0) 	   {return Struct(i)->laserOffsetW;}
  float       triggerTimeOffset(int i = 0)     {return 1e-6 * (mLaser ? laserOffset(i)  : offset(i));} // usec
  float       triggerTimeOffsetWest(int i = 0) {return 1e-6 * (mLaser ? laserOffsetW(i) :         0);} // usec
  void          SetLaser(bool k = true)         {mLaser = k;}
 private:
  bool        mLaser;
};
#endif
