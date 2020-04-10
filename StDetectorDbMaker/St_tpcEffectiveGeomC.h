#ifndef St_tpcEffectiveGeomC_h
#define St_tpcEffectiveGeomC_h

#include "tpcrs/config_structs.h"
#include "tpcEffectiveGeom.h"

struct St_tpcEffectiveGeomC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcEffectiveGeomC, tpcEffectiveGeom_st>
{
  Double_t 	drift_length_correction(Int_t i = 0) {return Struct(i)->drift_length_correction;}
  Double_t 	z_inner_offset(Int_t i = 0) 	  {return Struct(i)->z_inner_offset;}
  Double_t 	z_outer_offset(Int_t i = 0) 	  {return Struct(i)->z_outer_offset;}
  Double_t 	z_inner_offset_West(Int_t i = 0)  {return Struct(i)->z_inner_offset_West;}
  Double_t 	z_outer_offset_West(Int_t i = 0)  {return Struct(i)->z_outer_offset_West;}
  /*  Double_t 	scale(Int_t i = 0)                {return Struct(i)->scale;} */
};
#endif
