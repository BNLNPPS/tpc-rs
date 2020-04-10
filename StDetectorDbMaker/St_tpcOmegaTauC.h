#ifndef St_tpcOmegaTauC_h
#define St_tpcOmegaTauC_h

#include "tpcrs/config_structs.h"
#include "tpcOmegaTau.h"

struct St_tpcOmegaTauC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcOmegaTauC, tpcOmegaTau_st> {
 public:
  Float_t 	tensorV1(Int_t i = 0) 	     {return Struct(i)->tensorV1;}
  Float_t 	tensorV2(Int_t i = 0) 	     {return Struct(i)->tensorV2;}
  Float_t 	getOmegaTauTensorV1()        {return tensorV1();}
  Float_t 	getOmegaTauTensorV2()        {return tensorV2();}
  UInt_t        distortionCorrectionsMode(Int_t i = 0)
                                             {return Struct(i)->distortionCorrectionsMode;}
};
#endif
