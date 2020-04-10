#ifndef St_tpcAnodeHVC_h
#define St_tpcAnodeHVC_h

#include "tpcrs/config_structs.h"
#include "tpcAnodeHV.h"

struct St_tpcAnodeHVC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAnodeHVC, tpcAnodeHV_st>
{
  UShort_t 	 sector(Int_t i = 0) 	const {return Struct(i)->sector;}
  UShort_t 	 socket(Int_t i = 0) 	const {return Struct(i)->socket;}
  Float_t 	 voltage(Int_t i = 0) 	const;
  Bool_t	 livePadrow(Int_t sector = 1, Int_t padrow = 1) const { return voltagePadrow(sector, padrow) > 500; }
  Float_t	 voltagePadrow(Int_t sector = 1, Int_t padrow = 1) const ; // sector=1..24 , padrow=1..100
  Bool_t         tripped(Int_t sector = 1, Int_t padrow = 1) const { return (voltagePadrow(sector, padrow) < -100); }
  static  void   sockets(Int_t sector, Int_t padrow, Int_t &e1, Int_t &e2, Float_t &f2);
};
#endif
