#ifndef St_tpcGlobalPositionC_h
#define St_tpcGlobalPositionC_h

#include "tpcrs/config_structs.h"
#include "tpcGlobalPosition.h"

struct St_tpcGlobalPositionC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGlobalPositionC, tpcGlobalPosition_st>
{
  float 	LocalxShift(int i = 0)       const {return Struct(i)->LocalxShift;}
  float 	LocalyShift(int i = 0)       const {return Struct(i)->LocalyShift;}
  float 	LocalzShift(int i = 0)       const {return Struct(i)->LocalzShift;}
  /*  float 	PhiXY(int i = 0)  	       const {return Struct(i)->PhiXY;}	   */
  float 	PhiXZ(int i = 0)  	       const {return Struct(i)->PhiXZ;}
  float 	PhiYZ(int i = 0)  	       const {return Struct(i)->PhiYZ;}
  /*  float 	XX(int i = 0)  	       const {return Struct(i)->XX;}
      float 	YY(int i = 0)  	       const {return Struct(i)->YY;}
      float 	ZZ(int i = 0)  	       const {return Struct(i)->ZZ;}	    */
  float 	PhiXY_geom(int i = 0)        const {return Struct(i)->PhiXY_geom;}
  float 	PhiXZ_geom(int i = 0)        const {return Struct(i)->PhiXZ_geom;}
  float 	PhiYZ_geom(int i = 0)        const {return Struct(i)->PhiYZ_geom;}
  /*  float 	XX_geom(int i = 0)  	       const {return Struct(i)->XX_geom;}
      float 	YY_geom(int i = 0)  	       const {return Struct(i)->YY_geom;}
      float 	ZZ_geom(int i = 0)  	       const {return Struct(i)->ZZ_geom;}   */
  double  	TpcCenterPositionX()           const {return LocalxShift();}
  double  	TpcCenterPositionY()           const {return LocalyShift();}
  double  	TpcCenterPositionZ()           const {return LocalzShift();}
  double  	TpcRotationAroundGlobalAxisX() const {return PhiYZ_geom();}
  double  	TpcRotationAroundGlobalAxisY() const {return PhiXZ_geom();}
  double  	TpcRotationAroundGlobalAxisZ() const {return PhiXY_geom();}
  double  	TpcEFieldRotationX()           const {return PhiYZ();} /* YTWIST */
  double  	TpcEFieldRotationY() 	       const {return PhiXZ();} /* XTWIST */
  double      XTWIST()                       const {return  1e3 * TpcEFieldRotationY();}
  double      YTWIST()                       const {return -1e3 * TpcEFieldRotationX();}
  /* double  	TpcEFieldRotationZ() 	       const {return PhiXY();}              */
  double      X0()                           const {return LocalxShift();}
  double      Y0()                           const {return LocalyShift();}
  double      Z0()                           const {return LocalzShift();}
  double      alpha()                        const {return PhiYZ_geom();}
  double      beta()                         const {return PhiXZ_geom();}
  double      gamma()                        const {return PhiXY_geom();}
};
#endif
