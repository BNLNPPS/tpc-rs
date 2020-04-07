#ifndef StTpcSurveyC_h
#define StTpcSurveyC_h

#include "StDetectorDbMaker/St_SurveyC.h"
#include "tpcrs/enums.h"

class StTpcInnerSectorPosition : public St_SurveyC {// Inner part of sector to Super Sector
 public:
  static StTpcInnerSectorPosition* 	instance();
  StTpcInnerSectorPosition(St_Survey *table=0) : St_SurveyC(table) {}
  virtual ~StTpcInnerSectorPosition() {fgInstance = 0;}
 private:
  static StTpcInnerSectorPosition* fgInstance;
};
class StTpcOuterSectorPosition : public St_SurveyC {// Outer part of sector to Super Sector
 public:
  static StTpcOuterSectorPosition* 	instance();
  StTpcOuterSectorPosition(St_Survey *table=0) : St_SurveyC(table) {}
  virtual ~StTpcOuterSectorPosition() {fgInstance = 0;}
 private:
  static StTpcOuterSectorPosition* fgInstance;
};

class StTpcSuperSectorPosition : public St_SurveyC {// Extra rotation for whole Super Sector to half Tpc
 public:
  static StTpcSuperSectorPosition* 	instance();
  StTpcSuperSectorPosition(St_Survey *table=0) : St_SurveyC(table) {}
  virtual ~StTpcSuperSectorPosition() {fgInstance = 0;}
 private:
  static StTpcSuperSectorPosition* fgInstance;
};

class StTpcHalfPosition : public St_SurveyC {// Extra rotation for half of Tpc  to Tpc
 public:
  static StTpcHalfPosition* 	instance();
  StTpcHalfPosition(St_Survey *table=0) : St_SurveyC(table) {}
  virtual ~StTpcHalfPosition() {fgInstance = 0;}
  const TGeoHMatrix  &GetEastMatrix() {return  GetMatrix(TPC::Half::first);}
  const TGeoHMatrix  &GetWestMatrix() {return  GetMatrix(TPC::Half::second);}
  const TGeoHMatrix  &GetEastMatrixR() {return  GetMatrixR(TPC::Half::first);}
  const TGeoHMatrix  &GetWestMatrixR() {return  GetMatrixR(TPC::Half::second);}
  static void Normalize(TGeoHMatrix &R) {}
 private:
  static StTpcHalfPosition* fgInstance;
};

class StTpcPosition : public St_SurveyC {// Global position of TPC in Magnet
 public:
  static StTpcPosition* 	instance();
  StTpcPosition(St_Survey *table=0) : St_SurveyC(table) {}
  virtual ~StTpcPosition() {fgInstance = 0;}
  const TGeoHMatrix  &GetMatrix() {return  St_SurveyC::GetMatrix(0);}
 private:
  static StTpcPosition* fgInstance;
};
#endif
