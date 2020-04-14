#ifndef St_itpcPadGainT0C_h
#define St_itpcPadGainT0C_h

#include "tpcrs/config_structs.h"
#include "itpcPadGainT0.h"

struct St_itpcPadGainT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_itpcPadGainT0C, itpcPadGainT0_st>
{
  int 	run(int i = 0) 	const {return Struct(i)->run;}
  float 	Gain(int sector, int row, int pad) const
  {
    return ((sector > 0 && sector <= 24) && (row > 0 && row <= 40) && (pad > 0 && pad <= 120)) ?
           Struct()->Gain[sector - 1][row - 1][pad - 1] : 0;
  }
  float 	  T0(int sector, int row, int pad) const
  {
    return ((sector > 0 && sector <= 24) && (row > 0 && row <= 40) && (pad > 0 && pad <= 120)) ?
           Struct()->T0[sector - 1][row - 1][pad - 1] : 0;
  }
  bool    livePadrow(int sector, int row)
  {
    for (int pad = 1; pad <= 120; pad++) if (Gain(sector, row, pad) > 0) return kTRUE;

    return kFALSE;
  }
};
#endif
