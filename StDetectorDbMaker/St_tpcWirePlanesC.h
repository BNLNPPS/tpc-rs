#ifndef St_tpcWirePlanesC_h
#define St_tpcWirePlanesC_h

#include "tpcrs/config_structs.h"
#include "tpcWirePlanes.h"

struct St_tpcWirePlanesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcWirePlanesC, tpcWirePlanes_st>
{
  double 	anodeWireRadius(int i = 0) 	            {return Struct(i)->anodeWireRadius;}
  double 	frischGridWireRadius(int i = 0) 	    {return Struct(i)->frischGridWireRadius;}
  double 	gatingGridWireRadius(int i = 0) 	    {return Struct(i)->gatingGridWireRadius;}
  double 	anodeWirePitch(int i = 0) 	            {return Struct(i)->anodeWirePitch;}
  double 	frischGridWirePitch(int i = 0) 	    {return Struct(i)->frischGridWirePitch;}
  double 	gatingGridWirePitch(int i = 0) 	    {return Struct(i)->gatingGridWirePitch;}
  double 	innerSectorAnodeWirePadSep(int i = 0)     {return Struct(i)->innerSectorAnodeWirePadSep;}
  double 	innerSectorFrischGridPadSep(int i = 0)    {return Struct(i)->innerSectorFrischGridPadSep;}
  double 	innerSectorGatingGridPadSep(int i = 0)    {return Struct(i)->innerSectorGatingGridPadSep;}
  double 	outerSectorAnodeWirePadSep(int i = 0)     {return Struct(i)->outerSectorAnodeWirePadSep;}
  double 	outerSectorFrischGridPadSep(int i = 0)    {return Struct(i)->outerSectorFrischGridPadSep;}
  double 	outerSectorGatingGridPadSep(int i = 0)    {return Struct(i)->outerSectorGatingGridPadSep;}
  int 	numInnerSectorAnodeWires(int i = 0) 	    {return Struct(i)->numInnerSectorAnodeWires;}
  int 	numInnerSectorFrischGridWires(int i = 0)  {return Struct(i)->numInnerSectorFrischGridWires;}
  int 	numInnerSectorGatingGridWires(int i = 0)  {return Struct(i)->numInnerSectorGatingGridWires;}
  double 	firstInnerSectorAnodeWire(int i = 0) 	    {return Struct(i)->firstInnerSectorAnodeWire;}
  double 	firstInnerSectorFrischGridWire(int i = 0) {return Struct(i)->firstInnerSectorFrischGridWire;}
  double 	firstInnerSectorGatingGridWire(int i = 0) {return Struct(i)->firstInnerSectorGatingGridWire;}
  double 	lastInnerSectorAnodeWire(int i = 0) 	    {return Struct(i)->lastInnerSectorAnodeWire;}
  int 	numOuterSectorAnodeWires(int i = 0) 	    {return Struct(i)->numOuterSectorAnodeWires;}
  int 	numOuterSectorFrischGridWires(int i = 0)  {return Struct(i)->numOuterSectorFrischGridWires;}
  int 	numOuterSectorGatingGridWires(int i = 0)  {return Struct(i)->numOuterSectorGatingGridWires;}
  double 	firstOuterSectorAnodeWire(int i = 0) 	    {return Struct(i)->firstOuterSectorAnodeWire;}
  double 	firstOuterSectorFrischGridWire(int i = 0) {return Struct(i)->firstOuterSectorFrischGridWire;}
  double 	firstOuterSectorGatingGridWire(int i = 0) {return Struct(i)->firstOuterSectorGatingGridWire;}
  double 	lastOuterSectorAnodeWire(int i = 0) 	    {return Struct(i)->lastOuterSectorAnodeWire;}

  double      gateWireRadius(int i = 0)  {return gatingGridWireRadius(i);}
  double  	frischGridPitch(int i = 0) {return frischGridWirePitch(i);}
  double  	gatePitch(int i = 0)       {return gatingGridWirePitch(i);}

  double  	innerSectorAnodeWirePadPlaneSeparation(int i = 0)  {return innerSectorAnodeWirePadSep(i);}
  double  	innerSectorFrischGridPadPlaneSeparation(int i = 0) {return innerSectorFrischGridPadSep(i);}
  double  	innerSectorGatingGridPadPlaneSeparation(int i = 0) {return innerSectorGatingGridPadSep(i);}
  double  	outerSectorAnodeWirePadPlaneSeparation(int i = 0)  {return outerSectorAnodeWirePadSep(i);}
  double  	outerSectorFrischGridPadPlaneSeparation(int i = 0) {return outerSectorFrischGridPadSep(i);}
  double  	outerSectorGatingGridPadPlaneSeparation(int i = 0) {return outerSectorGatingGridPadSep(i);}

  int         numberOfInnerSectorAnodeWires(int i = 0)      {return numInnerSectorAnodeWires(i);}
  int   	numberOfInnerSectorFrischGridWires(int i = 0) {return numInnerSectorFrischGridWires(i);}
  int   	numberOfInnerSectorGatingGridWires(int i = 0) {return numInnerSectorGatingGridWires(i);}
  int   	numberOfOuterSectorAnodeWires(int i = 0)      {return numOuterSectorAnodeWires(i);}
  int   	numberOfOuterSectorFrischGridWires(int i = 0) {return numOuterSectorFrischGridWires(i);}
  int   	numberOfOuterSectorGatingGridWires(int i = 0) {return numOuterSectorGatingGridWires(i);}
};
#endif
