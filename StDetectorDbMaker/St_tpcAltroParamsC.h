#ifndef St_tpcAltroParamsC_h
#define St_tpcAltroParamsC_h

#include "tpcrs/config_structs.h"
#include "tpcAltroParams.h"

struct St_tpcAltroParamsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAltroParamsC, tpcAltroParams_st>
{
  Int_t 	N(Int_t i = 0) 	{return Struct(i)->N;}
  Int_t  	Threshold(Int_t i = 0) 	{return Struct(i)->Altro_thr;}
  Int_t  	MinSamplesaboveThreshold(Int_t i = 0) 	{return Struct(i)->Altro_seq;}
  Int_t  	K1(Int_t i = 0) 	{return Struct(i)->Altro_K1;}
  Int_t  	K2(Int_t i = 0) 	{return Struct(i)->Altro_K2;}
  Int_t  	K3(Int_t i = 0) 	{return Struct(i)->Altro_K3;}
  Int_t  	L1(Int_t i = 0) 	{return Struct(i)->Altro_L1;}
  Int_t  	L2(Int_t i = 0) 	{return Struct(i)->Altro_L2;}
  Int_t  	L3(Int_t i = 0) 	{return Struct(i)->Altro_L3;}
};
#endif
