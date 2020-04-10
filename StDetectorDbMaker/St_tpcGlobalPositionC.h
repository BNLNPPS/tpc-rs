#ifndef St_tpcGlobalPositionC_h
#define St_tpcGlobalPositionC_h

#include "tpcrs/config_structs.h"
#include "tpcGlobalPosition.h"

struct St_tpcGlobalPositionC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGlobalPositionC, tpcGlobalPosition_st>
{
  Float_t 	LocalxShift(Int_t i = 0)       const {return Struct(i)->LocalxShift;}
  Float_t 	LocalyShift(Int_t i = 0)       const {return Struct(i)->LocalyShift;}
  Float_t 	LocalzShift(Int_t i = 0)       const {return Struct(i)->LocalzShift;}
  /*  Float_t 	PhiXY(Int_t i = 0)  	       const {return Struct(i)->PhiXY;}	   */
  Float_t 	PhiXZ(Int_t i = 0)  	       const {return Struct(i)->PhiXZ;}
  Float_t 	PhiYZ(Int_t i = 0)  	       const {return Struct(i)->PhiYZ;}
  /*  Float_t 	XX(Int_t i = 0)  	       const {return Struct(i)->XX;}
      Float_t 	YY(Int_t i = 0)  	       const {return Struct(i)->YY;}
      Float_t 	ZZ(Int_t i = 0)  	       const {return Struct(i)->ZZ;}	    */
  Float_t 	PhiXY_geom(Int_t i = 0)        const {return Struct(i)->PhiXY_geom;}
  Float_t 	PhiXZ_geom(Int_t i = 0)        const {return Struct(i)->PhiXZ_geom;}
  Float_t 	PhiYZ_geom(Int_t i = 0)        const {return Struct(i)->PhiYZ_geom;}
  /*  Float_t 	XX_geom(Int_t i = 0)  	       const {return Struct(i)->XX_geom;}
      Float_t 	YY_geom(Int_t i = 0)  	       const {return Struct(i)->YY_geom;}
      Float_t 	ZZ_geom(Int_t i = 0)  	       const {return Struct(i)->ZZ_geom;}   */
  Double_t  	TpcCenterPositionX()           const {return LocalxShift();}
  Double_t  	TpcCenterPositionY()           const {return LocalyShift();}
  Double_t  	TpcCenterPositionZ()           const {return LocalzShift();}
  Double_t  	TpcRotationAroundGlobalAxisX() const {return PhiYZ_geom();}
  Double_t  	TpcRotationAroundGlobalAxisY() const {return PhiXZ_geom();}
  Double_t  	TpcRotationAroundGlobalAxisZ() const {return PhiXY_geom();}
  Double_t  	TpcEFieldRotationX()           const {return PhiYZ();} /* YTWIST */
  Double_t  	TpcEFieldRotationY() 	       const {return PhiXZ();} /* XTWIST */
  Double_t      XTWIST()                       const {return  1e3 * TpcEFieldRotationY();}
  Double_t      YTWIST()                       const {return -1e3 * TpcEFieldRotationX();}
  /* Double_t  	TpcEFieldRotationZ() 	       const {return PhiXY();}              */
  Double_t      X0()                           const {return LocalxShift();}
  Double_t      Y0()                           const {return LocalyShift();}
  Double_t      Z0()                           const {return LocalzShift();}
  Double_t      alpha()                        const {return PhiYZ_geom();}
  Double_t      beta()                         const {return PhiXZ_geom();}
  Double_t      gamma()                        const {return PhiXY_geom();}
};
#endif
