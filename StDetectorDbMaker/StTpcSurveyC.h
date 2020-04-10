#ifndef StTpcSurveyC_h
#define StTpcSurveyC_h

#include "StDetectorDbMaker/St_SurveyC.h"
#include "tpcrs/enums.h"
struct StTpcInnerSectorPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcInnerSectorPosition, Survey_st> {// Inner part of sector to Super Sector
};
struct StTpcOuterSectorPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcOuterSectorPosition, Survey_st> {// Outer part of sector to Super Sector
};

struct StTpcSuperSectorPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcSuperSectorPosition, Survey_st> {// Extra rotation for whole Super Sector to half Tpc
};

struct StTpcHalfPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcHalfPosition, Survey_st> {// Extra rotation for half of Tpc  to Tpc
  const TGeoHMatrix  &GetEastMatrix() {return  GetMatrix(TPC::Half::first);}
  const TGeoHMatrix  &GetWestMatrix() {return  GetMatrix(TPC::Half::second);}
  const TGeoHMatrix  &GetEastMatrixR() {return  GetMatrixR(TPC::Half::first);}
  const TGeoHMatrix  &GetWestMatrixR() {return  GetMatrixR(TPC::Half::second);}
  static void Normalize(TGeoHMatrix &R) {}
};

struct StTpcPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcPosition, Survey_st> {// Global position of TPC in Magnet
  const TGeoHMatrix  &GetMatrix() {return  St_SurveyC::GetMatrix(0);}
};
#endif
