#ifndef St_starMagOnlC_h
#define St_starMagOnlC_h

#include "tpcrs/config_structs.h"
#include "TMath.h"
#include "starMagOnl.h"

enum StMagnetPolarity {eUnknownMField, eFullMFieldPolB, eHalfMFieldPolB,
                       eZeroMField, eHalfMFieldPolA, eFullMFieldPolA
                      };

struct St_starMagOnlC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_starMagOnlC, starMagOnl_st>
{
  UInt_t 	runNumber(Int_t i = 0) 	{return Struct(i)->runNumber;}
  UInt_t 	time(Int_t i = 0) 	{return Struct(i)->time;}
  Double_t 	current(Int_t i = 0) 	{return Struct(i)->current;}
  Double_t      getScaleFactor(UInt_t time = 0) {return currentToScaleFactor(getMagnetCurrent(time));}
  Double_t      getMagnetCurrent(UInt_t time = 0)
  {
    if (! instance()) return 0;

    if (GetNRows() == 1 || time == 0) return current();

    Double_t tempCurrent = -9999;

    for (UInt_t i = 0; i < GetNRows() - 1; i++)
      if ( time >= getTimeEntry(i) && time <= getTimeEntry(i + 1) )
        if ( TMath::Abs(getMagnetCurrentEntry(i) - getMagnetCurrentEntry(i + 1)) < 50 )
          tempCurrent = getMagnetCurrentEntry(i);

    return tempCurrent;
  }
  StMagnetPolarity           getMagneticField(UInt_t time = 0)
  {
    StMagnetPolarity value = eUnknownMField;

    if (! instance()) return value;

    Double_t scaleFactor = getScaleFactor(time);

    if (scaleFactor == 1.0)	value = eFullMFieldPolA;

    if (scaleFactor == 0.5)	value = eHalfMFieldPolA;

    if (scaleFactor == 0.0)	value = eZeroMField;

    if (scaleFactor == -0.5)	value = eHalfMFieldPolB;

    if (scaleFactor == -1.0)	value = eFullMFieldPolB;

    return value;
  }
  UInt_t        getRunNumber() {return runNumber();}
  UInt_t        getTimeEntry(UInt_t i = 0) {return time(i);}
  Double_t      getMagnetCurrentEntry(UInt_t i = 0) {return current(i);}
  static Double_t  currentToScaleFactor(Double_t current)
  {
    Double_t value = -9999;

    if (! instance()) return value;

    if     (current < -4450 && current > -4550)	value = -1.0;
    else if (current < -2200 && current > -2300)	value = -0.5;
    else if (current >   -50 && current <    50)	value =  0.0;
    else if (current >  2200 && current <  2300)	value =  0.5;
    else if (current >  4450 && current <  4550)	value =  1.0;

    return value;
  }
};
#endif
