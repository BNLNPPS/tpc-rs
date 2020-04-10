#ifndef St_iTPCSurveyC_h
#define St_iTPCSurveyC_h

#include "tpcrs/config_structs.h"
#include "iTPCSurvey.h"

struct St_iTPCSurveyC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_iTPCSurveyC, iTPCSurvey_st> {
  Int_t 	Id(Int_t i = 0) 	const {return Struct(i)->Id;}
  Float_t 	Angle(Int_t i = 0) 	const {return Struct(i)->Angle;}
  Float_t 	dx(Int_t i = 0) 	const {return Struct(i)->dx;}
  Float_t 	dy(Int_t i = 0) 	const {return Struct(i)->dy;}
  Float_t 	ScaleX(Int_t i = 0) 	const {return Struct(i)->ScaleX;}
  Float_t 	ScaleY(Int_t i = 0) 	const {return Struct(i)->ScaleY;}
  Char_t* 	comment(Int_t i = 0) 	const {return Struct(i)->comment;}
};
#endif
