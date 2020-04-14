#ifndef St_iTPCSurveyC_h
#define St_iTPCSurveyC_h

#include "tpcrs/config_structs.h"
#include "iTPCSurvey.h"

struct St_iTPCSurveyC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_iTPCSurveyC, iTPCSurvey_st> {
  int 	Id(int i = 0) 	const {return Struct(i)->Id;}
  float 	Angle(int i = 0) 	const {return Struct(i)->Angle;}
  float 	dx(int i = 0) 	const {return Struct(i)->dx;}
  float 	dy(int i = 0) 	const {return Struct(i)->dy;}
  float 	ScaleX(int i = 0) 	const {return Struct(i)->ScaleX;}
  float 	ScaleY(int i = 0) 	const {return Struct(i)->ScaleY;}
  char* 	comment(int i = 0) 	const {return Struct(i)->comment;}
};
#endif
