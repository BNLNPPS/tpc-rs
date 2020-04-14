#ifndef St_tpcAnodeHVC_h
#define St_tpcAnodeHVC_h

#include "tpcrs/config_structs.h"
#include "tpcAnodeHV.h"

struct St_tpcAnodeHVC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAnodeHVC, tpcAnodeHV_st>
{
  unsigned short 	 sector(int i = 0) 	const {return Struct(i)->sector;}
  unsigned short 	 socket(int i = 0) 	const {return Struct(i)->socket;}
  float 	 voltage(int i = 0) 	const;
  bool	 livePadrow(int sector = 1, int padrow = 1) const { return voltagePadrow(sector, padrow) > 500; }
  float	 voltagePadrow(int sector = 1, int padrow = 1) const ; // sector=1..24 , padrow=1..100
  bool         tripped(int sector = 1, int padrow = 1) const { return (voltagePadrow(sector, padrow) < -100); }
  static  void   sockets(int sector, int padrow, int &e1, int &e2, float &f2);
};
#endif
