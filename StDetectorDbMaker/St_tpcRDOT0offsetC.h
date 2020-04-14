#ifndef St_tpcRDOT0offsetC_h
#define St_tpcRDOT0offsetC_h

#include "tpcrs/config_structs.h"
#include "tpcRDOT0offset.h"

struct St_tpcRDOT0offsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOT0offsetC, tpcRDOT0offset_st> {
  unsigned char* 	isShifted(int i = 0) 	const {return Struct(i)->isShifted;}
  bool        IsShfited(int sector) const {return isShifted()[sector-1];}
  float       T0(int sector, int padrow, int pad) const;
};
#endif
