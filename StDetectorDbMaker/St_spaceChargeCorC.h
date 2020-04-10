#ifndef St_spaceChargeCorC_h
#define St_spaceChargeCorC_h

#include "tpcrs/config_structs.h"
#include "spaceChargeCor.h"
#include "StDetectorDbMaker/StDetectorDbMagnet.h"

struct St_spaceChargeCorC : tpcrs::IConfigStruct {
  virtual spaceChargeCor_st* Struct(int i = 0) const = 0;
  Double_t 	fullFieldB(Int_t i = 0) 	{return Struct(i)->fullFieldB;}
  Double_t 	halfFieldB(Int_t i = 0) 	{return Struct(i)->halfFieldB;}
  Double_t 	zeroField(Int_t i = 0) 	        {return Struct(i)->zeroField;}
  Double_t 	halfFieldA(Int_t i = 0) 	{return Struct(i)->halfFieldA;}
  Double_t 	fullFieldA(Int_t i = 0) 	{return Struct(i)->fullFieldA;}
  Double_t 	satRate(Int_t i = 0) 	        {return Struct(i)->satRate;}
  Float_t 	factor(Int_t i = 0) 	        {return Struct(i)->factor;}
  Float_t 	detector(Int_t i = 0) 	        {return Struct(i)->detector;}
  Float_t 	offset(Int_t i = 0) 	        {return Struct(i)->offset;}
  Float_t 	getEWRatio(Int_t i = 0)	        {return Struct(i)->ewratio;}
  Double_t      getSpaceChargeCorrection(Double_t scaleFactor, Int_t i = 0){
    Double_t value = 0;
    if(scaleFactor < -.75 && scaleFactor > -1.25) value = fullFieldB(i);
    else if(scaleFactor < -0.25)	          value = halfFieldB(i);
    else if(scaleFactor < .25)	                  value = zeroField(i);
    else if(scaleFactor < 0.75)	                  value = halfFieldA(i);
    else if(scaleFactor < 1.25)	                  value = fullFieldA(i);
    return value;
  }
  Double_t getSpaceChargeCorrection(){return  getSpaceChargeCorrection(StDetectorDbMagnet::instance()->getScaleFactor());}
  Double_t getSpaceChargeCoulombs(Double_t scaleFactor);
  Double_t getSpaceChargeCoulombs(){return getSpaceChargeCoulombs(StDetectorDbMagnet::instance()->getScaleFactor());}
  Double_t getSpaceChargeSatRate(Int_t i = 0) {return satRate(i);}
  Float_t  getSpaceChargeFactor(Int_t i = 0)  {return factor(i);}
  Float_t  getSpaceChargeDetector(Int_t i = 0){return detector(i);}
  Float_t  getSpaceChargeOffset(Int_t i = 0)  {return offset(i);}

};

struct St_spaceChargeCorR1C : tpcrs::ConfigStruct<St_spaceChargeCorC, St_spaceChargeCorR1C, spaceChargeCor_st> {
};
struct St_spaceChargeCorR2C : tpcrs::ConfigStruct<St_spaceChargeCorC, St_spaceChargeCorR2C, spaceChargeCor_st> {
};
#endif
