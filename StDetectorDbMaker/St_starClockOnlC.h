#ifndef St_starClockOnlC_h
#define St_starClockOnlC_h

#include "tpcrs/config_structs.h"
#include "starClockOnl.h"

struct St_starClockOnlC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_starClockOnlC, starClockOnl_st>
{
  starClockOnl_st* 	Struct(Int_t i = 0);
  UInt_t    RunNumber(Int_t i = 0)           {return Struct(i)->runNumber;}
  Double_t  CurrentFrequency(Int_t i = 0)    {return Struct(i)->frequency;}
  UInt_t    Time(Int_t i = 0)                {return Struct(i)->time;}
  Double_t  Frequency(Int_t i = 0)           {return CurrentFrequency(i);}
  // depreciated
  UInt_t    getRunNumber(Int_t i = 0)        {return RunNumber(i);}
  Double_t  getCurrentFrequency(Int_t i = 0) {return CurrentFrequency(i);}
  UInt_t    getTime(Int_t i = 0)             {return Time(i);}
  Double_t  getFrequency(Int_t i = 0)        {return Frequency(i);}
  Double_t  samplingFrequency(Int_t i = 0)   {return 1e-6 * CurrentFrequency(i);}
};
#endif
