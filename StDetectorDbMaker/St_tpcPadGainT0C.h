#ifndef St_tpcPadGainT0C_h
#define St_tpcPadGainT0C_h

#include "tpcrs/config_structs.h"
#include "tpcPadGainT0.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
struct St_tpcPadGainT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadGainT0C, tpcPadGainT0_st>
{
  int 	run()           	const {return Struct()->run;}
  float 	Gain(int sector, int row, int pad) const
  {
    float gain = 0;

    if ((sector > 0 && sector <= 24) && (row > 0 && row <= St_tpcPadConfigC::instance()->padRows(sector)) && (pad > 0 && pad <= 182)) {
      gain = Struct()->Gain[sector - 1][row - 1][pad - 1];
    }

    return gain;
  }
  float 	  T0(int sector, int row, int pad) const
  {
    float t0 = 0;

    if ((sector > 0 && sector <= 24) && (row > 0 && row <= St_tpcPadConfigC::instance()->padRows(sector)) && (pad > 0 && pad <= 182)) {
      t0 = Struct()->T0[sector - 1][row - 1][pad - 1];
    }

    return t0;
  }
  bool    livePadrow(int sector, int row)
  {
    for (int pad = 1; pad <= 182; pad++) if (Gain(sector, row, pad) > 0) return kTRUE;

    return kFALSE;
  }
};
#endif
