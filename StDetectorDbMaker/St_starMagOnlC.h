#ifndef St_starMagOnlC_h
#define St_starMagOnlC_h

#include "tpcrs/config_structs.h"
#include "starMagOnl.h"

enum StMagnetPolarity {eUnknownMField, eFullMFieldPolB, eHalfMFieldPolB,
                       eZeroMField, eHalfMFieldPolA, eFullMFieldPolA
                      };

struct St_starMagOnlC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_starMagOnlC, starMagOnl_st>
{
  unsigned int 	runNumber(int i = 0) 	{return Struct(i)->runNumber;}
  unsigned int 	time(int i = 0) 	{return Struct(i)->time;}
  double 	current(int i = 0) 	{return Struct(i)->current;}
  double      getScaleFactor(unsigned int time = 0) {return currentToScaleFactor(getMagnetCurrent(time));}
  double      getMagnetCurrent(unsigned int time = 0)
  {
    if (! instance()) return 0;

    if (GetNRows() == 1 || time == 0) return current();

    double tempCurrent = -9999;

    for (unsigned int i = 0; i < GetNRows() - 1; i++)
      if ( time >= getTimeEntry(i) && time <= getTimeEntry(i + 1) )
        if ( std::abs(getMagnetCurrentEntry(i) - getMagnetCurrentEntry(i + 1)) < 50 )
          tempCurrent = getMagnetCurrentEntry(i);

    return tempCurrent;
  }
  StMagnetPolarity           getMagneticField(unsigned int time = 0)
  {
    StMagnetPolarity value = eUnknownMField;

    if (! instance()) return value;

    double scaleFactor = getScaleFactor(time);

    if (scaleFactor == 1.0)	value = eFullMFieldPolA;

    if (scaleFactor == 0.5)	value = eHalfMFieldPolA;

    if (scaleFactor == 0.0)	value = eZeroMField;

    if (scaleFactor == -0.5)	value = eHalfMFieldPolB;

    if (scaleFactor == -1.0)	value = eFullMFieldPolB;

    return value;
  }
  unsigned int        getRunNumber() {return runNumber();}
  unsigned int        getTimeEntry(unsigned int i = 0) {return time(i);}
  double      getMagnetCurrentEntry(unsigned int i = 0) {return current(i);}
  static double  currentToScaleFactor(double current)
  {
    double value = -9999;

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
