#ifndef St_tpcAnodeHVavgC_h
#define St_tpcAnodeHVavgC_h

#include "tpcrs/config_structs.h"
#include "tpcAnodeHVavg.h"

struct St_tpcAnodeHVavgC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAnodeHVavgC, tpcAnodeHVavg_st>
{
  unsigned short          sector(int i = 0) 	const {return Struct(i)->sector;}
  unsigned short          socket(int i = 0) 	const {return Struct(i)->socket;}
  float 	    voltage(int i = 0) 	const;
  float 	    rms(int i = 0) 	        const {return Struct(i)->rms;}
  int 	    numentries(int i = 0) 	const {return Struct(i)->numentries;}
  int 	    numoutliers(int i = 0) 	const {return Struct(i)->numoutliers;}
  bool	    livePadrow(int sec = 1, int padrow = 1) const { return voltagePadrow(sec, padrow) > 500; }
  float	    voltagePadrow(int sec = 1, int padrow = 1) const; // sector=1..24 , padrow=1..100
  bool            tripped(int sec = 1, int padrow = 1)       const;// { return (voltage() < -100); }
};
#endif
