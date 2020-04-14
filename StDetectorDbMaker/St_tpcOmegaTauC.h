#ifndef St_tpcOmegaTauC_h
#define St_tpcOmegaTauC_h

#include "tpcrs/config_structs.h"
#include "tpcOmegaTau.h"

struct St_tpcOmegaTauC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcOmegaTauC, tpcOmegaTau_st> {
 public:
  float 	tensorV1(int i = 0) 	     {return Struct(i)->tensorV1;}
  float 	tensorV2(int i = 0) 	     {return Struct(i)->tensorV2;}
  float 	getOmegaTauTensorV1()        {return tensorV1();}
  float 	getOmegaTauTensorV2()        {return tensorV2();}
  unsigned int        distortionCorrectionsMode(int i = 0)
                                             {return Struct(i)->distortionCorrectionsMode;}
};
#endif
