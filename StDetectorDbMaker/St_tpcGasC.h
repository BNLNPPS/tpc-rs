#ifndef St_tpcGasC_h
#define St_tpcGasC_h

#include "tpcrs/config_structs.h"
#include "tpcGas.h"

struct St_tpcGasC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGasC, tpcGas_st>
{
  float 	barometricPressure(int i = 0) {return Struct(i)->barometricPressure;}
  float 	inputTPCGasPressure(int i = 0) {return Struct(i)->inputTPCGasPressure;}
  float 	nitrogenPressure(int i = 0) 	{return Struct(i)->nitrogenPressure;}
  float 	gasPressureDiff(int i = 0) 	{return Struct(i)->gasPressureDiff;}
  float 	inputGasTemperature(int i = 0) {return Struct(i)->inputGasTemperature;}
  float 	outputGasTemperature(int i = 0) {return Struct(i)->outputGasTemperature;}
  float 	flowRateArgon1(int i = 0) 	{return Struct(i)->flowRateArgon1;}
  float 	flowRateArgon2(int i = 0) 	{return Struct(i)->flowRateArgon2;}
  float 	flowRateMethane(int i = 0) 	{return Struct(i)->flowRateMethane;}
  float 	percentMethaneIn(int i = 0) 	{return Struct(i)->percentMethaneIn;}
  float 	ppmOxygenIn(int i = 0) 	{return Struct(i)->ppmOxygenIn;}
  float 	flowRateExhaust(int i = 0) 	{return Struct(i)->flowRateExhaust;}
  float 	percentMethaneOut(int i = 0) 	{return Struct(i)->percentMethaneOut;}
  float 	ppmWaterOut(int i = 0) 	{return Struct(i)->ppmWaterOut;}
  float 	ppmOxygenOut(int i = 0) 	{return Struct(i)->ppmOxygenOut;}
  float 	flowRateRecirculation(int i = 0) {return Struct(i)->flowRateRecirculation;}
};
#endif
