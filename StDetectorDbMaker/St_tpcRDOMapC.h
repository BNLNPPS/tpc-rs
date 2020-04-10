#ifndef St_tpcRDOMapC_h
#define St_tpcRDOMapC_h

#include "tpcrs/config_structs.h"
#include "tpcRDOMap.h"

struct St_tpcRDOMapC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOMapC, tpcRDOMap_st> {
  UChar_t 	nrows(Int_t i = 0) 	const {return Struct(i)->nrows;}
  UChar_t 	index(Int_t i = 0) 	const {return Struct(i)->idx;}
  UChar_t 	row(Int_t i = 0) 	const {return Struct(i)->row;}
  UChar_t 	padMin(Int_t i = 0) 	const {return Struct(i)->padMin;}
  UChar_t 	padMax(Int_t i = 0) 	const {return Struct(i)->padMax;}
  UChar_t 	rdoI(Int_t i = 0) 	const {return Struct(i)->rdo;}
  Int_t         rdo(Int_t padrow, Int_t pad = 0) const;
};
#endif
