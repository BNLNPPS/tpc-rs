#ifndef St_tpcAnodeHVavgC_h
#define St_tpcAnodeHVavgC_h

#include "tpcrs/config_structs.h"
#include "tpcAnodeHVavg.h"

struct St_tpcAnodeHVavgC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAnodeHVavgC, tpcAnodeHVavg_st>
{
  UShort_t          sector(Int_t i = 0) 	const {return Struct(i)->sector;}
  UShort_t          socket(Int_t i = 0) 	const {return Struct(i)->socket;}
  Float_t 	    voltage(Int_t i = 0) 	const;
  Float_t 	    rms(Int_t i = 0) 	        const {return Struct(i)->rms;}
  Int_t 	    numentries(Int_t i = 0) 	const {return Struct(i)->numentries;}
  Int_t 	    numoutliers(Int_t i = 0) 	const {return Struct(i)->numoutliers;}
  Bool_t	    livePadrow(Int_t sec = 1, Int_t padrow = 1) const { return voltagePadrow(sec, padrow) > 500; }
  Float_t	    voltagePadrow(Int_t sec = 1, Int_t padrow = 1) const; // sector=1..24 , padrow=1..100
  Bool_t            tripped(Int_t sec = 1, Int_t padrow = 1)       const;// { return (voltage() < -100); }
};
#endif
