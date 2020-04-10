#ifndef St_tpcPadResponseC_h
#define St_tpcPadResponseC_h

#include "tpcrs/config_structs.h"
#include "tpcPadResponse.h"

struct St_tpcPadResponseC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadResponseC, tpcPadResponse_st>
{
  Float_t 	innerGasGainFluctuation(Int_t i = 0) 	{return Struct(i)->innerGasGainFluctuation;}
  Float_t 	outerGasGainFluctuation(Int_t i = 0) 	{return Struct(i)->outerGasGainFluctuation;}
  Float_t 	innerPadResponseSigma(Int_t i = 0) 	{return Struct(i)->innerPadResponseSigma;}
  Float_t 	outerPadResponseSigma(Int_t i = 0) 	{return Struct(i)->outerPadResponseSigma;}
  Float_t 	innerWirePadCoupling(Int_t i = 0) 	{return Struct(i)->innerWirePadCoupling;}
  Float_t 	outerWirePadCoupling(Int_t i = 0) 	{return Struct(i)->outerWirePadCoupling;}
  Float_t 	innerRowNormalization(Int_t i = 0) 	{return Struct(i)->innerRowNormalization;}
  Float_t 	outerRowNormalization(Int_t i = 0) 	{return Struct(i)->outerRowNormalization;}
  Float_t* 	BoundaryOfStepFunctions(Int_t i = 0) 	{return Struct(i)->BoundaryOfStepFunctions;}
  Float_t* 	innerChargeFractionConstants(Int_t i = 0) {return Struct(i)->innerChargeFractionConstants;}
  Float_t* 	outerChargeFractionConstants(Int_t i = 0) {return Struct(i)->outerChargeFractionConstants;}
  Float_t 	errorFunctionRange(Int_t i = 0) 	{return Struct(i)->errorFunctionRange;}
  Int_t 	errorFunctionEntry(Int_t i = 0) 	{return Struct(i)->errorFunctionEntry;}
  Float_t 	longitudinalDiffusionConstant(Int_t i = 0) {return Struct(i)->longitudinalDiffusionConstant;}
  Float_t 	transverseDiffusionConstant(Int_t i = 0) {return Struct(i)->transverseDiffusionConstant;}
  Float_t 	InnerOuterFactor(Int_t i = 0) 	        {return Struct(i)->InnerOuterFactor;}
};
#endif
