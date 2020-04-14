#ifndef St_TpcAvgCurrentC_h
#define St_TpcAvgCurrentC_h

#include "tpcrs/config_structs.h"
#include "TpcAvgCurrent.h"
#include "StDetectorDbMaker/St_tpcAnodeHVC.h"
struct St_TpcAvgCurrentC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcAvgCurrentC, TpcAvgCurrent_st>
{
  int 	run(int i = 0) 	const {return Struct(i)->run;}
  int 	start_time(int i = 0) 	const {return Struct(i)->start_time;}
  int 	stop_time(int i = 0) 	const {return Struct(i)->stop_time;}
  static int  ChannelFromRow(int sector, int row);
  static int  ChannelFromSocket(int socket);
  float       AvCurrent(int sector = 1, int channel = 1);
  /* {
     return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
     Struct()->AvCurrent[8*(sector-1)+channel-1] :
     0;} */
  float       AvCurrSocket(int sector = 1, int socket = 1) {return AvCurrent(sector, ChannelFromSocket(socket));}
  float       AvCurrRow(int sector = 1, int row = 1) {return AvCurrent(sector, ChannelFromRow(sector, row));}
  float       AcCharge(int sector = 1, int channel = 1);
  /* {
     return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
     Struct()->AcCharge[8*(sector-1)+channel-1] :
     0;
     } */
  float       AcChargeSocket(int sector = 1, int socket = 1) {return AcCharge(sector, ChannelFromSocket(socket));}
  float       AcChargeRow(int sector = 1, int row = 1) {return AcCharge(sector, ChannelFromRow(sector, row));}
  float       AcChargeL(int sector = 1, int channel = 1); // C/cm
  float       AcChargeRowL(int sector = 1, int row = 1) {return AcChargeL(sector, ChannelFromRow(sector, row));}
};
#endif
