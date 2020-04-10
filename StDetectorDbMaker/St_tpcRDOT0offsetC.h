#ifndef St_tpcRDOT0offsetC_h
#define St_tpcRDOT0offsetC_h

#include "tpcrs/config_structs.h"
#include "tpcRDOT0offset.h"

struct St_tpcRDOT0offsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOT0offsetC, tpcRDOT0offset_st> {
  UChar_t* 	isShifted(Int_t i = 0) 	const {return Struct(i)->isShifted;}
  Bool_t        IsShfited(Int_t sector) const {return isShifted()[sector-1];}
  Float_t       T0(Int_t sector, Int_t padrow, Int_t pad) const;
};
#endif
