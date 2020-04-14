#ifndef St_tpcGridLeakC_h
#define St_tpcGridLeakC_h

#include "tpcrs/config_structs.h"
#include "tpcGridLeak.h"

enum StGLpos {
  kGLinner=0,
  kGLmiddl=1,
  kGLouter=2
};
struct St_tpcGridLeakC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGridLeakC, tpcGridLeak_st> {
  double 	InnerGLRadius(int i = 0) 	{return Struct(i)->InnerGLRadius;}
  double 	MiddlGLRadius(int i = 0) 	{return Struct(i)->MiddlGLRadius;}
  double 	OuterGLRadius(int i = 0) 	{return Struct(i)->OuterGLRadius;}
  double 	InnerGLWidth(int i = 0) 	{return Struct(i)->InnerGLWidth;}
  double 	MiddlGLWidth(int i = 0) 	{return Struct(i)->MiddlGLWidth;}
  double 	OuterGLWidth(int i = 0) 	{return Struct(i)->OuterGLWidth;}
  double 	InnerGLStrength(int i = 0) 	{return Struct(i)->InnerGLStrength;}
  double 	MiddlGLStrength(int i = 0) 	{return Struct(i)->MiddlGLStrength;}
  double 	OuterGLStrength(int i = 0) 	{return Struct(i)->OuterGLStrength;}
  double      getGridLeakStrength(StGLpos pos){
    switch (pos) {
      case (kGLinner) : return InnerGLStrength();
      case (kGLmiddl) : return MiddlGLStrength();
      case (kGLouter) : return OuterGLStrength();
    }
    return 0;
  }
  double      getGridLeakRadius(StGLpos pos) {
    switch (pos) {
      case (kGLinner) : return InnerGLRadius();
      case (kGLmiddl) : return MiddlGLRadius();
      case (kGLouter) : return OuterGLRadius();
    }
    return 0;
  }
  double      getGridLeakWidth(StGLpos pos) {
    switch (pos) {
      case (kGLinner) : return InnerGLWidth();
      case (kGLmiddl) : return MiddlGLWidth();
      case (kGLouter) : return OuterGLWidth();
    }
    return 0;
  }
};
#endif
