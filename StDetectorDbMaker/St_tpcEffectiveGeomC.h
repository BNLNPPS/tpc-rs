#ifndef St_tpcEffectiveGeomC_h
#define St_tpcEffectiveGeomC_h

#include "tpcrs/config_structs.h"
#include "tpcEffectiveGeom.h"

struct St_tpcEffectiveGeomC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcEffectiveGeomC, tpcEffectiveGeom_st>
{
  double 	drift_length_correction(int i = 0) {return Struct(i)->drift_length_correction;}
  double 	z_inner_offset(int i = 0) 	  {return Struct(i)->z_inner_offset;}
  double 	z_outer_offset(int i = 0) 	  {return Struct(i)->z_outer_offset;}
  double 	z_inner_offset_West(int i = 0)  {return Struct(i)->z_inner_offset_West;}
  double 	z_outer_offset_West(int i = 0)  {return Struct(i)->z_outer_offset_West;}
  /*  double 	scale(int i = 0)                {return Struct(i)->scale;} */
};
#endif
