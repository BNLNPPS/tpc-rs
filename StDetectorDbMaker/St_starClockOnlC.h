#ifndef St_starClockOnlC_h
#define St_starClockOnlC_h

#include "tpcrs/config_structs.h"
#include "starClockOnl.h"

struct St_starClockOnlC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_starClockOnlC, starClockOnl_st>
{
  unsigned int    RunNumber(int i = 0)           {return Struct(i)->runNumber;}
  double  CurrentFrequency(int i = 0)    {return Struct(i)->frequency;}
  unsigned int    Time(int i = 0)                {return Struct(i)->time;}
  double  Frequency(int i = 0)           {return CurrentFrequency(i);}
  // depreciated
  unsigned int    getRunNumber(int i = 0)        {return RunNumber(i);}
  double  getCurrentFrequency(int i = 0) {return CurrentFrequency(i);}
  unsigned int    getTime(int i = 0)             {return Time(i);}
  double  getFrequency(int i = 0)        {return Frequency(i);}
  double  samplingFrequency(int i = 0)   {return 1e-6 * CurrentFrequency(i);}
};
#endif
