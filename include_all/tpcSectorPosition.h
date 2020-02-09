/* tpcSectorPosition.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcSectorPosition.idl

  Table: tpcSectorPosition

       description:
      This table rotates the inner and outer sector points using the following
      steps:
      1) move x position (sector 12 coordinates) by outerSectorLocalxShift
      2) Rotate outer points clockwise by outerSectorRotationAngle about
         point (0,123cm).
      3) Rotate inner points clockwise by innerSectorRotationAngle about
         point (0,123cm)

 */
#ifndef TPCSECTORPOSITION_H
#define TPCSECTORPOSITION_H
#define TPCSECTORPOSITION_SPEC \
"struct tpcSectorPosition { \
	float innerSectorLocalxShift; \
	float innerSectorLocalyShift; \
	float innerSectorRotationAngle; \
	float innerSectorCovMatrix; \
	float outerSectorLocalxShift; \
	float outerSectorLocalyShift; \
	float outerSectorRotationAngle; \
	float outerSectorCovMatrix; \
};"
typedef struct tpcSectorPosition_st {
	float innerSectorLocalxShift; /*   cm : shift in local x coord.  */
	float innerSectorLocalyShift; /*   cm : shift in local y coord.  */
	float innerSectorRotationAngle; /*   degrees : clockwise rotation  */
	float innerSectorCovMatrix; /*   0  */
	float outerSectorLocalxShift; /*   cm : shift in local x coord.  */
	float outerSectorLocalyShift; /*   cm : shift in local y coord.  */
	float outerSectorRotationAngle; /*   degrees : clockwise rotation  */
	float outerSectorCovMatrix; /*   0  */
} TPCSECTORPOSITION_ST;
#endif /* TPCSECTORPOSITION_H */
