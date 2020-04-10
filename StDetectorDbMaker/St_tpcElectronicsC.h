#ifndef St_tpcElectronicsC_h
#define St_tpcElectronicsC_h

#include "tpcrs/config_structs.h"
#include "tpcElectronics.h"
#include "StDetectorDbMaker/St_starClockOnlC.h"
struct St_tpcElectronicsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcElectronicsC, tpcElectronics_st>
{
  Int_t 	numberOfTimeBins(Int_t i = 0) 	{return Struct(i)->numberOfTimeBins;}
  Double_t 	nominalGain(Int_t i = 0) 	{return Struct(i)->nominalGain;}
  //  Double_t 	samplingFrequency(Int_t i = 0) 	{return Struct(i)->samplingFrequency;}  obsolete
  Double_t      samplingFrequency(Int_t i = 0) {return 1e-6 * St_starClockOnlC::instance()->CurrentFrequency(i);}
  Double_t 	tZero(Int_t i = 0) 	        {return Struct(i)->tZero;}
  Double_t 	adcCharge(Int_t i = 0) 	        {return Struct(i)->adcCharge;}
  Double_t 	adcConversion(Int_t i = 0) 	{return Struct(i)->adcConversion;}
  Double_t 	averagePedestal(Int_t i = 0) 	{return Struct(i)->averagePedestal;}
  Double_t 	shapingTime(Int_t i = 0) 	{return Struct(i)->shapingTime;}
  Double_t 	tau(Int_t i = 0) 	        {return Struct(i)->tau;}
};
#endif
