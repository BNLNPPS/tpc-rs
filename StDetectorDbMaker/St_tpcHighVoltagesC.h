#ifndef St_tpcHighVoltagesC_h
#define St_tpcHighVoltagesC_h

#include "tpcrs/config_structs.h"
#include "tpcHighVoltages.h"

struct St_tpcHighVoltagesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcHighVoltagesC, tpcHighVoltages_st> {
  Float_t 	cathode(Int_t i = 0)          {return Struct(i)->cathode;}
  Float_t 	gatedGridRef(Int_t i = 0)     {return Struct(i)->gatedGridRef;}
  Float_t* 	gridLeakWallTip(Int_t i = 0)  {return Struct(i)->gridLeakWallTip;}
  Float_t* 	gridLeakWallSide(Int_t i = 0) {return Struct(i)->gridLeakWallSide;}
  Double_t      getCathodeVoltage()           {return cathode();}
  Double_t      getGGVoltage()                {return gatedGridRef();}
  Double_t      getGridLeakWallTip(Int_t sector = 1)  {return gridLeakWallTip()[sector-1];}
  Double_t      getGridLeakWallSide(Int_t sector = 1) {return gridLeakWallSide()[sector-1];}
};
#endif
