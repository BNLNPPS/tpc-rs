#ifndef St_tpcAltroParamsC_h
#define St_tpcAltroParamsC_h

#include "tpcrs/config_structs.h"
#include "tpcAltroParams.h"

struct St_tpcAltroParamsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAltroParamsC, tpcAltroParams_st>
{
  int 	N(int i = 0) 	{return Struct(i)->N;}
  int  	Threshold(int i = 0) 	{return Struct(i)->Altro_thr;}
  int  	MinSamplesaboveThreshold(int i = 0) 	{return Struct(i)->Altro_seq;}
  int  	K1(int i = 0) 	{return Struct(i)->Altro_K1;}
  int  	K2(int i = 0) 	{return Struct(i)->Altro_K2;}
  int  	K3(int i = 0) 	{return Struct(i)->Altro_K3;}
  int  	L1(int i = 0) 	{return Struct(i)->Altro_L1;}
  int  	L2(int i = 0) 	{return Struct(i)->Altro_L2;}
  int  	L3(int i = 0) 	{return Struct(i)->Altro_L3;}
};
#endif
