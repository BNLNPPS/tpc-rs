#ifndef St_tpcHighVoltagesC_h
#define St_tpcHighVoltagesC_h

#include "tpcrs/config_structs.h"
#include "tpcHighVoltages.h"

struct St_tpcHighVoltagesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcHighVoltagesC, tpcHighVoltages_st> {
  float 	cathode(int i = 0)          {return Struct(i)->cathode;}
  float 	gatedGridRef(int i = 0)     {return Struct(i)->gatedGridRef;}
  float* 	gridLeakWallTip(int i = 0)  {return Struct(i)->gridLeakWallTip;}
  float* 	gridLeakWallSide(int i = 0) {return Struct(i)->gridLeakWallSide;}
  double      getCathodeVoltage()           {return cathode();}
  double      getGGVoltage()                {return gatedGridRef();}
  double      getGridLeakWallTip(int sector = 1)  {return gridLeakWallTip()[sector-1];}
  double      getGridLeakWallSide(int sector = 1) {return gridLeakWallSide()[sector-1];}
};
#endif
