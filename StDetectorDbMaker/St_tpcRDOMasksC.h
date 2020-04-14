#ifndef St_tpcRDOMasksC_h
#define St_tpcRDOMasksC_h

#include "tpcrs/config_structs.h"
#include "tpcRDOMasks.h"
#include "StDetectorDbMaker/St_tpcPadPlanesC.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
struct St_tpcRDOMasksC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOMasksC, tpcRDOMasks_st>
{
  unsigned int 	runNumber(int i = 0) 	        {return Struct(i)->runNumber;}
  unsigned int 	sector(int i = 0) 	        {return Struct(i)->sector;}
  unsigned int 	mask(int i = 0) 	        {return Struct(i)->mask;}
  unsigned int        getSectorMask(unsigned int sector);
  static unsigned int rdoForPadrow(int row)   //Function returns the rdo board number for a given padrow index. Range of map used is 1-45.
  {
    unsigned int rdo = 0;

    if      (row > 0 && row <=  8) rdo = 1;
    else if (row > 8 && row <= 13) rdo = 2;
    else if (row > 13 && row <= 21) rdo = 3;
    else if (row > 21 && row <= 29) rdo = 4;
    else if (row > 29 && row <= 37) rdo = 5;
    else if (row > 37 && row <= 45) rdo = 6;

    return rdo;
  }
  static unsigned int rdoForPadrow(int sector, int row)   //Function returns the rdo board number for a given padrow index. Range of map used is 1-45.
  {
    if (St_tpcPadConfigC::instance()->iTpc(sector)) return 8;

    return rdoForPadrow(row);
  }
  bool        isOn(int sector, int rdo)
  {
    if (St_tpcPadConfigC::instance()->iTpc(sector)) return 1;

    if (sector < 1 || sector > 24 || rdo < 1 || rdo > 6)	return 0;

    unsigned int MASK = getSectorMask(sector);
    MASK = MASK >> (rdo - 1);
    MASK &= 0x00000001;
    return MASK;
  }
  bool       isRowOn(int sector, int row) {return isOn(sector, rdoForPadrow(sector, row));}
};
#endif
