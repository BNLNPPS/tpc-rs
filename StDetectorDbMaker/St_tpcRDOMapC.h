#ifndef St_tpcRDOMapC_h
#define St_tpcRDOMapC_h

#include "tpcrs/config_structs.h"
#include "tpcRDOMap.h"

struct St_tpcRDOMapC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOMapC, tpcRDOMap_st> {
  unsigned char 	nrows(int i = 0) 	const {return Struct(i)->nrows;}
  unsigned char 	index(int i = 0) 	const {return Struct(i)->idx;}
  unsigned char 	row(int i = 0) 	const {return Struct(i)->row;}
  unsigned char 	padMin(int i = 0) 	const {return Struct(i)->padMin;}
  unsigned char 	padMax(int i = 0) 	const {return Struct(i)->padMax;}
  unsigned char 	rdoI(int i = 0) 	const {return Struct(i)->rdo;}
  int         rdo(int padrow, int pad = 0) const;
};
#endif
