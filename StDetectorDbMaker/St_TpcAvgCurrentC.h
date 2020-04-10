#ifndef St_TpcAvgCurrentC_h
#define St_TpcAvgCurrentC_h

#include "tpcrs/config_structs.h"
#include "TpcAvgCurrent.h"
#include "StDetectorDbMaker/St_tpcAnodeHVC.h"
struct St_TpcAvgCurrentC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcAvgCurrentC, TpcAvgCurrent_st>
{
  Int_t 	run(Int_t i = 0) 	const {return Struct(i)->run;}
  Int_t 	start_time(Int_t i = 0) 	const {return Struct(i)->start_time;}
  Int_t 	stop_time(Int_t i = 0) 	const {return Struct(i)->stop_time;}
  static Int_t  ChannelFromRow(Int_t sector, Int_t row);
  static Int_t  ChannelFromSocket(Int_t socket);
  Float_t       AvCurrent(Int_t sector = 1, Int_t channel = 1);
  /* {
     return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
     Struct()->AvCurrent[8*(sector-1)+channel-1] :
     0;} */
  Float_t       AvCurrSocket(Int_t sector = 1, Int_t socket = 1) {return AvCurrent(sector, ChannelFromSocket(socket));}
  Float_t       AvCurrRow(Int_t sector = 1, Int_t row = 1) {return AvCurrent(sector, ChannelFromRow(sector, row));}
  Float_t       AcCharge(Int_t sector = 1, Int_t channel = 1);
  /* {
     return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
     Struct()->AcCharge[8*(sector-1)+channel-1] :
     0;
     } */
  Float_t       AcChargeSocket(Int_t sector = 1, Int_t socket = 1) {return AcCharge(sector, ChannelFromSocket(socket));}
  Float_t       AcChargeRow(Int_t sector = 1, Int_t row = 1) {return AcCharge(sector, ChannelFromRow(sector, row));}
  Float_t       AcChargeL(Int_t sector = 1, Int_t channel = 1); // C/cm
  Float_t       AcChargeRowL(Int_t sector = 1, Int_t row = 1) {return AcChargeL(sector, ChannelFromRow(sector, row));}
};
#endif
