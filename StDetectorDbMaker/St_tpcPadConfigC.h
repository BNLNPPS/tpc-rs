#ifndef St_tpcPadConfigC_h
#define St_tpcPadConfigC_h

#include "tpcrs/config_structs.h"
#include "tpcPadConfig.h"

struct St_tpcPadConfigC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadConfigC, tpcPadConfig_st>
{
  unsigned char          iTpc(int sector);
  unsigned char          iTPC(int sector) {return iTpc(sector);}
  int 	   padRows(int sector);
  int 	   innerPadRows(int sector);
  int 	   innerPadRows48(int sector);
  int 	   innerPadRows52(int sector);
  int 	   outerPadRows(int sector);
  int 	   superInnerPadRows(int sector);
  int 	   superOuterPadRows(int sector);
  double 	   innerSectorPadWidth(int sector);
  double 	   innerSectorPadLength(int sector);
  double 	   innerSectorPadPitch(int sector);
  double 	   innerSectorRowPitch1(int sector);
  double 	   innerSectorRowPitch2(int sector);
  double 	   firstPadRow(int sector);
  double 	   firstOuterSectorPadRow(int sector);
  double 	   lastOuterSectorPadRow(int sector);
  double 	   firstRowWidth(int sector);
  double 	   lastRowWidth(int sector);
  double 	   outerSectorPadWidth(int sector);
  double 	   outerSectorPadLength(int sector);
  double 	   outerSectorPadPitch(int sector);
  double 	   outerSectorRowPitch(int sector);
  double 	   outerSectorLength(int sector);
  double 	   ioSectorSeparation(int sector);
  double 	   innerSectorEdge(int sector);
  double 	   outerSectorEdge(int sector);
  double 	   innerSectorPadPlaneZ(int sector);
  double 	   outerSectorPadPlaneZ(int sector);
  int* 	   innerPadsPerRow(int sector);
  int* 	   outerPadsPerRow(int sector);
  int            padsPerRow(int sector, int row = 1);
  double* 	   innerRowRadii(int sector);
  double* 	   outerRowRadii(int sector);
  //               taken from StRItpcPadPlane
  int            numberOfRows(int sector);
  int            numberOfInnerRows(int sector);
  int            numberOfInnerRows48(int sector);
  int            numberOfInnerRows52(int sector);
  int            numberOfOuterRows(int sector);
  bool           isRowInRange(int sector, int row);
  double         radialDistanceAtRow(int sector, int row);
  int            numberOfPadsAtRow(int sector, int row);
  double         PadWidthAtRow(int sector, int row);
  double 	   PadLengthAtRow(int sector, int row);
  double 	   PadPitchAtRow(int sector, int row);
  double 	   RowPitchAtRow(int sector, int row);
  int            indexForRowPad(int sector, int row, int pad);
  bool             isiTpcSector(int sector) { return iTpc(sector) == 1; }
  bool             isiTpcPadRow(int sector, int row) { return iTpc(sector) && row >= 1 && row <= numberOfInnerRows(sector); }
  bool             isInnerPadRow(int sector, int row) { return row <= numberOfInnerRows(sector); }
  int            IsRowInner(int sector, int row) {return (row <= innerPadRows(sector)) ? 1 : 0;}
};
#endif
