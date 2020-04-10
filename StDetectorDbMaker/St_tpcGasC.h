#ifndef St_tpcGasC_h
#define St_tpcGasC_h

#include "tpcrs/config_structs.h"
#include "tpcGas.h"

struct St_tpcGasC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGasC, tpcGas_st>
{
  Float_t 	barometricPressure(Int_t i = 0) {return Struct(i)->barometricPressure;}
  Float_t 	inputTPCGasPressure(Int_t i = 0) {return Struct(i)->inputTPCGasPressure;}
  Float_t 	nitrogenPressure(Int_t i = 0) 	{return Struct(i)->nitrogenPressure;}
  Float_t 	gasPressureDiff(Int_t i = 0) 	{return Struct(i)->gasPressureDiff;}
  Float_t 	inputGasTemperature(Int_t i = 0) {return Struct(i)->inputGasTemperature;}
  Float_t 	outputGasTemperature(Int_t i = 0) {return Struct(i)->outputGasTemperature;}
  Float_t 	flowRateArgon1(Int_t i = 0) 	{return Struct(i)->flowRateArgon1;}
  Float_t 	flowRateArgon2(Int_t i = 0) 	{return Struct(i)->flowRateArgon2;}
  Float_t 	flowRateMethane(Int_t i = 0) 	{return Struct(i)->flowRateMethane;}
  Float_t 	percentMethaneIn(Int_t i = 0) 	{return Struct(i)->percentMethaneIn;}
  Float_t 	ppmOxygenIn(Int_t i = 0) 	{return Struct(i)->ppmOxygenIn;}
  Float_t 	flowRateExhaust(Int_t i = 0) 	{return Struct(i)->flowRateExhaust;}
  Float_t 	percentMethaneOut(Int_t i = 0) 	{return Struct(i)->percentMethaneOut;}
  Float_t 	ppmWaterOut(Int_t i = 0) 	{return Struct(i)->ppmWaterOut;}
  Float_t 	ppmOxygenOut(Int_t i = 0) 	{return Struct(i)->ppmOxygenOut;}
  Float_t 	flowRateRecirculation(Int_t i = 0) {return Struct(i)->flowRateRecirculation;}
};
#endif
