#ifndef St_trgTimeOffsetC_h
#define St_trgTimeOffsetC_h

#include "tpcrs/config_structs.h"
#include "trgTimeOffset.h"

struct St_trgTimeOffsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_trgTimeOffsetC, trgTimeOffset_st>
{
  Float_t 	offset(Int_t i = 0)     	   {return Struct(i)->offset;}
  Float_t 	laserOffset(Int_t i = 0) 	   {return Struct(i)->laserOffset;}
  Float_t 	laserOffsetW(Int_t i = 0) 	   {return Struct(i)->laserOffsetW;}
  Float_t       triggerTimeOffset(Int_t i = 0)     {return 1e-6 * (mLaser ? laserOffset(i)  : offset(i));} // usec
  Float_t       triggerTimeOffsetWest(Int_t i = 0) {return 1e-6 * (mLaser ? laserOffsetW(i) :         0);} // usec
  void          SetLaser(Bool_t k = kTRUE)         {mLaser = k;}
 private:
  Bool_t        mLaser;
};
#endif
