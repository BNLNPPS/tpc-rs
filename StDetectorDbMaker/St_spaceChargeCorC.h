#ifndef St_spaceChargeCorC_h
#define St_spaceChargeCorC_h

#include "tpcrs/config_structs.h"
#include "spaceChargeCor.h"
#include "StDetectorDbMaker/StDetectorDbMagnet.h"

struct St_spaceChargeCorC : tpcrs::IConfigStruct {
  virtual spaceChargeCor_st* Struct(int i = 0) const = 0;
  double 	fullFieldB(int i = 0) 	{return Struct(i)->fullFieldB;}
  double 	halfFieldB(int i = 0) 	{return Struct(i)->halfFieldB;}
  double 	zeroField(int i = 0) 	        {return Struct(i)->zeroField;}
  double 	halfFieldA(int i = 0) 	{return Struct(i)->halfFieldA;}
  double 	fullFieldA(int i = 0) 	{return Struct(i)->fullFieldA;}
  double 	satRate(int i = 0) 	        {return Struct(i)->satRate;}
  float 	factor(int i = 0) 	        {return Struct(i)->factor;}
  float 	detector(int i = 0) 	        {return Struct(i)->detector;}
  float 	offset(int i = 0) 	        {return Struct(i)->offset;}
  float 	getEWRatio(int i = 0)	        {return Struct(i)->ewratio;}
  double      getSpaceChargeCorrection(double scaleFactor, int i = 0){
    double value = 0;
    if(scaleFactor < -.75 && scaleFactor > -1.25) value = fullFieldB(i);
    else if(scaleFactor < -0.25)	          value = halfFieldB(i);
    else if(scaleFactor < .25)	                  value = zeroField(i);
    else if(scaleFactor < 0.75)	                  value = halfFieldA(i);
    else if(scaleFactor < 1.25)	                  value = fullFieldA(i);
    return value;
  }
  double getSpaceChargeCorrection(){return  getSpaceChargeCorrection(StDetectorDbMagnet::instance()->getScaleFactor());}
  double getSpaceChargeCoulombs(double scaleFactor);
  double getSpaceChargeCoulombs(){return getSpaceChargeCoulombs(StDetectorDbMagnet::instance()->getScaleFactor());}
  double getSpaceChargeSatRate(int i = 0) {return satRate(i);}
  float  getSpaceChargeFactor(int i = 0)  {return factor(i);}
  float  getSpaceChargeDetector(int i = 0){return detector(i);}
  float  getSpaceChargeOffset(int i = 0)  {return offset(i);}

};

struct St_spaceChargeCorR1C : tpcrs::ConfigStruct<St_spaceChargeCorC, St_spaceChargeCorR1C, spaceChargeCor_st> {
};
struct St_spaceChargeCorR2C : tpcrs::ConfigStruct<St_spaceChargeCorC, St_spaceChargeCorR2C, spaceChargeCor_st> {
};
#endif
