#ifndef St_tpcElectronicsC_h
#define St_tpcElectronicsC_h

#include "tpcrs/config_structs.h"
#include "tpcElectronics.h"
#include "StDetectorDbMaker/St_starClockOnlC.h"
struct St_tpcElectronicsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcElectronicsC, tpcElectronics_st>
{
  int 	numberOfTimeBins(int i = 0) 	{return Struct(i)->numberOfTimeBins;}
  double 	nominalGain(int i = 0) 	{return Struct(i)->nominalGain;}
  //  double 	samplingFrequency(int i = 0) 	{return Struct(i)->samplingFrequency;}  obsolete
  double      samplingFrequency(int i = 0) {return 1e-6 * St_starClockOnlC::instance()->CurrentFrequency(i);}
  double 	tZero(int i = 0) 	        {return Struct(i)->tZero;}
  double 	adcCharge(int i = 0) 	        {return Struct(i)->adcCharge;}
  double 	adcConversion(int i = 0) 	{return Struct(i)->adcConversion;}
  double 	averagePedestal(int i = 0) 	{return Struct(i)->averagePedestal;}
  double 	shapingTime(int i = 0) 	{return Struct(i)->shapingTime;}
  double 	tau(int i = 0) 	        {return Struct(i)->tau;}
};
#endif
