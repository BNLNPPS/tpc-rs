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
  Double_t 	InnerGLRadius(Int_t i = 0) 	{return Struct(i)->InnerGLRadius;}
  Double_t 	MiddlGLRadius(Int_t i = 0) 	{return Struct(i)->MiddlGLRadius;}
  Double_t 	OuterGLRadius(Int_t i = 0) 	{return Struct(i)->OuterGLRadius;}
  Double_t 	InnerGLWidth(Int_t i = 0) 	{return Struct(i)->InnerGLWidth;}
  Double_t 	MiddlGLWidth(Int_t i = 0) 	{return Struct(i)->MiddlGLWidth;}
  Double_t 	OuterGLWidth(Int_t i = 0) 	{return Struct(i)->OuterGLWidth;}
  Double_t 	InnerGLStrength(Int_t i = 0) 	{return Struct(i)->InnerGLStrength;}
  Double_t 	MiddlGLStrength(Int_t i = 0) 	{return Struct(i)->MiddlGLStrength;}
  Double_t 	OuterGLStrength(Int_t i = 0) 	{return Struct(i)->OuterGLStrength;}
  Double_t      getGridLeakStrength(StGLpos pos){
    switch (pos) {
      case (kGLinner) : return InnerGLStrength();
      case (kGLmiddl) : return MiddlGLStrength();
      case (kGLouter) : return OuterGLStrength();
    }
    return 0;
  }
  Double_t      getGridLeakRadius(StGLpos pos) {
    switch (pos) {
      case (kGLinner) : return InnerGLRadius();
      case (kGLmiddl) : return MiddlGLRadius();
      case (kGLouter) : return OuterGLRadius();
    }
    return 0;
  }
  Double_t      getGridLeakWidth(StGLpos pos) {
    switch (pos) {
      case (kGLinner) : return InnerGLWidth();
      case (kGLmiddl) : return MiddlGLWidth();
      case (kGLouter) : return OuterGLWidth();
    }
    return 0;
  }
};
#endif
