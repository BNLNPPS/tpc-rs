#ifndef St_tpcPadResponseC_h
#define St_tpcPadResponseC_h

#include "tpcrs/config_structs.h"
#include "tpcPadResponse.h"

struct St_tpcPadResponseC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadResponseC, tpcPadResponse_st>
{
  float 	innerGasGainFluctuation(int i = 0) 	{return Struct(i)->innerGasGainFluctuation;}
  float 	outerGasGainFluctuation(int i = 0) 	{return Struct(i)->outerGasGainFluctuation;}
  float 	innerPadResponseSigma(int i = 0) 	{return Struct(i)->innerPadResponseSigma;}
  float 	outerPadResponseSigma(int i = 0) 	{return Struct(i)->outerPadResponseSigma;}
  float 	innerWirePadCoupling(int i = 0) 	{return Struct(i)->innerWirePadCoupling;}
  float 	outerWirePadCoupling(int i = 0) 	{return Struct(i)->outerWirePadCoupling;}
  float 	innerRowNormalization(int i = 0) 	{return Struct(i)->innerRowNormalization;}
  float 	outerRowNormalization(int i = 0) 	{return Struct(i)->outerRowNormalization;}
  float* 	BoundaryOfStepFunctions(int i = 0) 	{return Struct(i)->BoundaryOfStepFunctions;}
  float* 	innerChargeFractionConstants(int i = 0) {return Struct(i)->innerChargeFractionConstants;}
  float* 	outerChargeFractionConstants(int i = 0) {return Struct(i)->outerChargeFractionConstants;}
  float 	errorFunctionRange(int i = 0) 	{return Struct(i)->errorFunctionRange;}
  int 	errorFunctionEntry(int i = 0) 	{return Struct(i)->errorFunctionEntry;}
  float 	longitudinalDiffusionConstant(int i = 0) {return Struct(i)->longitudinalDiffusionConstant;}
  float 	transverseDiffusionConstant(int i = 0) {return Struct(i)->transverseDiffusionConstant;}
  float 	InnerOuterFactor(int i = 0) 	        {return Struct(i)->InnerOuterFactor;}
};
#endif
