#ifndef St_tpcPadPlanesC_h
#define St_tpcPadPlanesC_h

#include "tpcrs/config_structs.h"
#include "tpcPadPlanes.h"

struct St_tpcPadPlanesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadPlanesC, tpcPadPlanes_st>
{
  int 	padRows(int i = 0) 	         {return Struct(i)->padRows;}
  int 	innerPadRows(int i = 0) 	 {return Struct(i)->innerPadRows;}
  int 	innerPadRows48(int i = 0) 	 {return Struct(i)->innerPadRows48;}
  int 	innerPadRows52(int i = 0) 	 {return Struct(i)->innerPadRows52;}
  int 	outerPadRows(int i = 0) 	 {return Struct(i)->outerPadRows;}
  int 	superInnerPadRows(int i = 0) 	 {return Struct(i)->superInnerPadRows;}
  int 	superOuterPadRows(int i = 0) 	 {return Struct(i)->superOuterPadRows;}
  double 	innerSectorPadWidth(int i = 0) {return Struct(i)->innerSectorPadWidth;}
  double 	innerSectorPadLength(int i = 0) {return Struct(i)->innerSectorPadLength;}
  double 	innerSectorPadPitch(int i = 0) {return Struct(i)->innerSectorPadPitch;}
  double 	innerSectorRowPitch1(int i = 0) {return Struct(i)->innerSectorRowPitch1;}
  double 	innerSectorRowPitch2(int i = 0) {return Struct(i)->innerSectorRowPitch2;}
  double 	firstPadRow(int i = 0) 	 {return Struct(i)->firstPadRow;}
  double 	firstOuterSectorPadRow(int i = 0) {return Struct(i)->firstOuterSectorPadRow;}
  double 	lastOuterSectorPadRow(int i = 0) {return Struct(i)->lastOuterSectorPadRow;}
  double 	firstRowWidth(int i = 0) 	 {return Struct(i)->firstRowWidth;}
  double 	lastRowWidth(int i = 0) 	 {return Struct(i)->lastRowWidth;}
  double 	outerSectorPadWidth(int i = 0) {return Struct(i)->outerSectorPadWidth;}
  double 	outerSectorPadLength(int i = 0) {return Struct(i)->outerSectorPadLength;}
  double 	outerSectorPadPitch(int i = 0) {return Struct(i)->outerSectorPadPitch;}
  double 	outerSectorRowPitch(int i = 0) {return Struct(i)->outerSectorRowPitch;}
  double 	outerSectorLength(int i = 0) 	 {return Struct(i)->outerSectorLength;}
  double 	ioSectorSeparation(int i = 0)  {return Struct(i)->ioSectorSeparation;}
  double 	innerSectorEdge(int i = 0) 	 {return Struct(i)->innerSectorEdge;}
  double 	outerSectorEdge(int i = 0) 	 {return Struct(i)->outerSectorEdge;}
  double 	innerSectorPadPlaneZ(int i = 0) {return Struct(i)->innerSectorPadPlaneZ;}
  double 	outerSectorPadPlaneZ(int i = 0) {return Struct(i)->outerSectorPadPlaneZ;}
  int* 	innerPadsPerRow(int i = 0) 	 {return Struct(i)->innerPadsPerRow;}
  int* 	outerPadsPerRow(int i = 0) 	 {return Struct(i)->outerPadsPerRow;}
  int         padsPerRow(int row = 1)
  {
    return (row <= innerPadRows()) ?
           innerPadsPerRow()[row - 1] :
           outerPadsPerRow()[row - 1 - innerPadRows()];
  }
  double* 	innerRowRadii(int i = 0) 	 {return Struct(i)->innerRowRadii;}
  double* 	outerRowRadii(int i = 0) 	 {return Struct(i)->outerRowRadii;}
  // taken from StRTpcPadPlane
  int         numberOfRows()                   {return padRows();}
  int         numberOfInnerRows()              {return innerPadRows();}
  int         numberOfInnerRows48()            {return innerPadRows48();}
  int         numberOfInnerRows52()            {return innerPadRows52();}
  int         numberOfOuterRows()              {return outerPadRows();}
  bool        isRowInRange(int row)          {return (row >= 1 && row <= numberOfRows()) ? kTRUE : kFALSE;}
  double      radialDistanceAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows() ) return innerRowRadii()[row - 1];
    else                            return outerRowRadii()[row - 1 - numberOfInnerRows()];
  }
  int   numberOfPadsAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows() ) return innerPadsPerRow()[row - 1];

    return outerPadsPerRow()[row - 1 - numberOfInnerRows()];
  }
  double PadWidthAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows()) return innerSectorPadWidth();

    return outerSectorPadWidth();
  }
  double PadLengthAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows()) return innerSectorPadLength();

    return outerSectorPadLength();
  }
  double PadPitchAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows()) return innerSectorPadPitch();

    return outerSectorPadPitch();
  }
  double RowPitchAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows48() ) return innerSectorRowPitch1();
    else if (row > numberOfInnerRows48() && row <= numberOfInnerRows()) return innerSectorRowPitch2();

    return outerSectorRowPitch();
  }
  int indexForRowPad(int row, int pad)
  {
    if (pad > numberOfPadsAtRow(row)) return -1;

    int index = 0;

    if (row > 0 && row <= numberOfInnerRows() )             for (int i = 1; i < row; i++) index += numberOfPadsAtRow(i);
    else if (row > numberOfInnerRows() && row <= numberOfRows()) for (int i = numberOfInnerRows() + 1; i < row; i++)  index += numberOfPadsAtRow(i);

    index += pad - 1;
    return index;
  }
};
#endif
