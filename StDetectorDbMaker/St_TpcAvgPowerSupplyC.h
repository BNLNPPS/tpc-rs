#ifndef St_TpcAvgPowerSupplyC_h
#define St_TpcAvgPowerSupplyC_h

#include "tpcrs/config_structs.h"
#include "TpcAvgPowerSupply.h"
#include "StDetectorDbMaker/St_TpcAvgCurrentC.h"
struct St_TpcAvgPowerSupplyC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcAvgPowerSupplyC, TpcAvgPowerSupply_st>
{
  int 	run(int i = 0) 	const {return Struct(i)->run;}
  int 	start_time(int i = 0) 	const {return Struct(i)->start_time;}
  int 	stop_time(int i = 0) 	const {return Struct(i)->stop_time;}
  float* 	Current(int i = 0) 	const {return Struct(i)->Current;}
  float* 	Charge(int i = 0) 	const {return Struct(i)->Charge;}
  float* 	Voltage(int i = 0) 	const {return Struct(i)->Voltage;}
  float	voltagePadrow(int sec = 1, int padrow = 1) const; // sector=1..24 , padrow=1..100
  bool        tripped(int sec = 1, int row = 1) const {return voltagePadrow(sec, row) < -100;}
  static int  ChannelFromRow(int sector, int row) {return St_TpcAvgCurrentC::ChannelFromRow(sector, row);}
  static int  ChannelFromSocket(int socket) {return St_TpcAvgCurrentC::ChannelFromSocket(socket);}
  float       AvCurrent(int sector = 1, int channel = 1)
  {
    return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
           Struct()->Current[8 * (sector - 1) + channel - 1] :
           0;
  }
  float       AvCurrSocket(int sector = 1, int socket = 1) {return AvCurrent(sector, ChannelFromSocket(socket));}
  float       AvCurrRow(int sector = 1, int row = 1) {return AvCurrent(sector, ChannelFromRow(sector, row));}
  float       AcCharge(int sector = 1, int channel = 1)
  {
    return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
           Struct()->Charge[8 * (sector - 1) + channel - 1] :
           0;
  }
  float       AcChargeSocket(int sector = 1, int socket = 1) {return AcCharge(sector, ChannelFromSocket(socket));}
  float       AcChargeRow(int sector = 1, int row = 1) {return AcCharge(sector, ChannelFromRow(sector, row));}
  float       AcChargeL(int sector = 1, int channel = 1); // C/cm
  float       AcChargeRowL(int sector = 1, int row = 1) {return AcChargeL(sector, ChannelFromRow(sector, row));}
  bool        livePadrow(int sec = 1, int padrow = 1) const { return voltagePadrow(sec, padrow) >  500;}
};
#endif
