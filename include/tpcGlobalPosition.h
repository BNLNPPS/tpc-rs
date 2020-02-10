/* tpcGlobalPosition.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcGlobalPosition.idl

  Table: tpcGlobalPosition

       description:

 */
#ifndef TPCGLOBALPOSITION_H
#define TPCGLOBALPOSITION_H
#define TPCGLOBALPOSITION_SPEC \
"struct tpcGlobalPosition { \
	float LocalxShift; \
	float LocalyShift; \
	float LocalzShift; \
	float PhiXY; \
	float PhiXZ; \
	float PhiYZ; \
	float XX; \
	float YY; \
	float ZZ; \
	float PhiXY_geom; \
	float PhiXZ_geom; \
	float PhiYZ_geom; \
	float XX_geom; \
	float YY_geom; \
	float ZZ_geom; \
};"
typedef struct tpcGlobalPosition_st {
	float LocalxShift; /* cm : x position of TPC center in magnet frame  */
	float LocalyShift; /* cm : y position of TPC center in magnet frame  */
	float LocalzShift; /* cm : z position of TPC center in magnet frame  */
	float PhiXY; /* radians: rotation angle around z axis  (not used) */
	float PhiXZ; /* radians: rotation angle around y axis  XTWIST */
	float PhiYZ; /* radians: rotation angle around x axis  YTWIST */
	float XX; /* XX element of rotation matrix  (not used) */
	float YY; /* YY element of rotation matrix  (not used) */
	float ZZ; /* ZZ element of rotation matrix  (not used) */
	float PhiXY_geom; /* radians: geometrical rotation angle around z axis psi,  -gamma  (not used) */
	float PhiXZ_geom; /* radians: geometrical rotation angle around y axis theta,-beta  */
	float PhiYZ_geom; /* radians: geometrical rotation angle around x axis psi,  -alpha */
	float XX_geom; /* XX element of geometrical rotation matrix  (not used) */
	float YY_geom; /* YY element of geometrical rotation matrix  (not used) */
	float ZZ_geom; /* ZZ element of geometrical rotation matrix  (not used) */
} TPCGLOBALPOSITION_ST;
#endif /* TPCGLOBALPOSITION_H */
