#ifndef St_itpcPadGainT0C_h
#define St_itpcPadGainT0C_h

#include "tpcrs/config_structs.h"
#include "itpcPadGainT0.h"

struct St_itpcPadGainT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_itpcPadGainT0C, itpcPadGainT0_st>
{
  Int_t 	run(Int_t i = 0) 	const {return Struct(i)->run;}
  Float_t 	Gain(Int_t sector, Int_t row, Int_t pad) const
  {
    return ((sector > 0 && sector <= 24) && (row > 0 && row <= 40) && (pad > 0 && pad <= 120)) ?
           Struct()->Gain[sector - 1][row - 1][pad - 1] : 0;
  }
  Float_t 	  T0(Int_t sector, Int_t row, Int_t pad) const
  {
    return ((sector > 0 && sector <= 24) && (row > 0 && row <= 40) && (pad > 0 && pad <= 120)) ?
           Struct()->T0[sector - 1][row - 1][pad - 1] : 0;
  }
  Bool_t    livePadrow(Int_t sector, Int_t row)
  {
    for (Int_t pad = 1; pad <= 120; pad++) if (Gain(sector, row, pad) > 0) return kTRUE;

    return kFALSE;
  }
};
#endif
