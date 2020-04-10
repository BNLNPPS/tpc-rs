#ifndef St_tpcSCGLC_h
#define St_tpcSCGLC_h
#include "tpcrs/config_structs.h"
#include "tpcSCGL.h"


struct St_tpcSCGLC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcSCGLC, tpcSCGL_st> {
  float* SC()                             {return Struct()->SC;}
  float* SCoffset()                       {return Struct()->SCoffset;}
  float* SCexponent()                     {return Struct()->SCexponent;}
  float* SCscaler()                       {return Struct()->SCscaler;}
  float* GL()                             {return Struct()->GL;}
  float* GLoffset()                       {return Struct()->GLoffset;}
  float  GLradius()                       {return Struct()->GLradius;}
  float  GLwidth()                        {return Struct()->GLwidth;}
  int    mode()                           {return Struct()->mode;}
  char*  comment()                        {return Struct()->comment;}
};
#endif
