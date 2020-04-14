#ifndef St_tpcDimensionsC_h
#define St_tpcDimensionsC_h

#include "tpcrs/config_structs.h"
#include "tpcDimensions.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
#include "StDetectorDbMaker/St_tpcEffectiveGeomC.h"
#include "StDetectorDbMaker/St_tpcWirePlanesC.h"
struct St_tpcDimensionsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcDimensionsC, tpcDimensions_st>
{
  int 	numberOfSectors(int i = 0) 	{return Struct(i)->numberOfSectors;}
  double 	tpcInnerRadius(int i = 0) 	{return Struct(i)->tpcInnerRadius;}
  double 	tpcOuterRadius(int i = 0) 	{return Struct(i)->tpcOuterRadius;}
  double 	tpcTotalLength(int i = 0) 	{return Struct(i)->tpcTotalLength;}
  double 	wheelInnerRadius(int i = 0) 	{return Struct(i)->wheelInnerRadius;}
  double 	wheelOuterRadius(int i = 0) 	{return Struct(i)->wheelOuterRadius;}
  double 	wheelThickness(int i = 0) 	{return Struct(i)->wheelThickness;}
  double 	senseGasOuterRadius(int i = 0) {return Struct(i)->senseGasOuterRadius;}
  double 	tpeaThickness(int i = 0) 	{return Struct(i)->tpeaThickness;}
  double 	cathodeInnerRadius(int i = 0) {return Struct(i)->cathodeInnerRadius;}
  double 	cathodeOuterRadius(int i = 0) {return Struct(i)->cathodeOuterRadius;}
  double 	cathodeThickness(int i = 0) 	{return Struct(i)->cathodeThickness;}
  double 	outerCuThickness(int i = 0) 	{return Struct(i)->outerCuThickness;}
  double 	outerKaptonThickness(int i = 0) {return Struct(i)->outerKaptonThickness;}
  double 	outerNomexThickness(int i = 0) {return Struct(i)->outerNomexThickness;}
  double 	outerGlueThickness(int i = 0) {return Struct(i)->outerGlueThickness;}
  double 	outerInsGasThickness(int i = 0) {return Struct(i)->outerInsGasThickness;}
  double 	outerAlThickness(int i = 0) 	{return Struct(i)->outerAlThickness;}
  double 	outerAlHoneycombThickness(int i = 0) 	{return Struct(i)->outerAlHoneycombThickness;}
  double 	innerGlueThickness(int i = 0) {return Struct(i)->innerGlueThickness;}
  double 	innerNomexThickness(int i = 0) {return Struct(i)->innerNomexThickness;}
  double 	innerKaptonThickness(int i = 0) {return Struct(i)->innerKaptonThickness;}
  double 	innerAlThickness(int i = 0) 	{return Struct(i)->innerAlThickness;}
  double 	innerGapWidI(int i = 0) 	{return Struct(i)->innerGapWidI;}
  double 	innerGapWidO(int i = 0) 	{return Struct(i)->innerGapWidO;}
  double 	innerGapHeit(int i = 0) 	{return Struct(i)->innerGapHeit;}
  double 	innerGapRad(int i = 0) 	{return Struct(i)->innerGapRad;}
  double 	innerInWidth(int i = 0) 	{return Struct(i)->innerInWidth;}
  double 	innerOutWidth(int i = 0) 	{return Struct(i)->innerOutWidth;}
  double 	innerHeight(int i = 0) 	{return Struct(i)->innerHeight;}
  double 	innerPPDepth(int i = 0) 	{return Struct(i)->innerPPDepth;}
  double 	innerAlDepth(int i = 0) 	{return Struct(i)->innerAlDepth;}
  double 	innerMWCDepth(int i = 0) 	{return Struct(i)->innerMWCDepth;}
  double 	innerBoundary(int i = 0) 	{return Struct(i)->innerBoundary;}
  double 	innerRCenter(int i = 0) 	{return Struct(i)->innerRCenter;}
  double 	innerMWCInn(int i = 0) 	{return Struct(i)->innerMWCInn;}
  double 	innerMWCOut(int i = 0) 	{return Struct(i)->innerMWCOut;}
  double 	innerMVCHei(int i = 0) 	{return Struct(i)->innerMVCHei;}
  int 	innerAirGaps(int i = 0) 	{return Struct(i)->innerAirGaps;}
  int 	innerExtraAl(int i = 0) 	{return Struct(i)->innerExtraAl;}
  double* 	innerZGaps(int i = 0) 	{return Struct(i)->innerZGaps;}
  double* 	innerZGapsSize(int i = 0) 	{return Struct(i)->innerZGapsSize;}
  double* 	innerXExtraAl(int i = 0) 	{return Struct(i)->innerXExtraAl;}
  double* 	innerZExtraAl(int i = 0) 	{return Struct(i)->innerZExtraAl;}
  double* 	innerDXExtraAl(int i = 0) 	{return Struct(i)->innerDXExtraAl;}
  double* 	innerDZExtraAl(int i = 0) 	{return Struct(i)->innerDZExtraAl;}
  double 	outerGapWidI(int i = 0) 	{return Struct(i)->outerGapWidI;}
  double 	outerGapWidO(int i = 0) 	{return Struct(i)->outerGapWidO;}
  double 	outerGapHeit(int i = 0) 	{return Struct(i)->outerGapHeit;}
  double 	outerGapRad(int i = 0) 	{return Struct(i)->outerGapRad;}
  double 	outerInWidth(int i = 0) 	{return Struct(i)->outerInWidth;}
  double 	outerOutWidth(int i = 0) 	{return Struct(i)->outerOutWidth;}
  double 	outerHeight(int i = 0) 	{return Struct(i)->outerHeight;}
  double 	outerPPDepth(int i = 0) 	{return Struct(i)->outerPPDepth;}
  double 	outerAlDepth(int i = 0) 	{return Struct(i)->outerAlDepth;}
  double 	outerMWCDepth(int i = 0) 	{return Struct(i)->outerMWCDepth;}
  double 	outerBoundary(int i = 0) 	{return Struct(i)->outerBoundary;}
  double 	outerRCenter(int i = 0) 	{return Struct(i)->outerRCenter;}
  double 	outerMWCInn(int i = 0) 	{return Struct(i)->outerMWCInn;}
  double 	outerMWCOut(int i = 0) 	{return Struct(i)->outerMWCOut;}
  double 	outerMVCHei(int i = 0) 	{return Struct(i)->outerMVCHei;}
  int 	outerAirGaps(int i = 0) 	{return Struct(i)->outerAirGaps;}
  int 	outerExtraAl(int i = 0) 	{return Struct(i)->outerExtraAl;}
  double* 	outerZGaps(int i = 0) 	{return Struct(i)->outerZGaps;}
  double* 	outerZGapsSize(int i = 0) 	{return Struct(i)->outerZGapsSize;}
  double* 	outerXExtraAl(int i = 0) 	{return Struct(i)->outerXExtraAl;}
  double* 	outerZExtraAl(int i = 0) 	{return Struct(i)->outerZExtraAl;}
  double* 	outerDXExtraAl(int i = 0) 	{return Struct(i)->outerDXExtraAl;}
  double* 	outerDZExtraAl(int i = 0) 	{return Struct(i)->outerDZExtraAl;}
  double      gatingGridZ(int sector = 20)
  {
    return St_tpcPadConfigC::instance()->outerSectorPadPlaneZ(sector)
           - St_tpcWirePlanesC::instance()->outerSectorGatingGridPadPlaneSeparation();
  }
  double      zInnerOffset()                  {return St_tpcEffectiveGeomC::instance()->z_inner_offset();}
  double      zOuterOffset()                  {return St_tpcEffectiveGeomC::instance()->z_outer_offset();}
  double      zInnerOffset_West()             {return St_tpcEffectiveGeomC::instance()->z_inner_offset_West();}
  double      zOuterOffset_West()             {return St_tpcEffectiveGeomC::instance()->z_outer_offset_West();}
  //TPC field cage parameters:
  double ifcRadius() {return tpcInnerRadius();}
  double ofcRadius() {return tpcOuterRadius();}
};
#endif
