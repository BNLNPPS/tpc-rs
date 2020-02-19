/* tpcSCGL.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcSCGL.idl

  Table: tpcSCGL

       description: Table for SpaceCharge and GridLeak Correction parameters
       SC parameters: 4 array elements for west, plus 4 for east  8,
          each of the 4 elements are additive terms in the SC formula
       GL parameters: 24 array elements for TPC sectors,
          plus radius and width of the charge sheet
       scaler and mode definitions in StDetectorDbMakerSt_tpcSCGLC.h

 */
#ifndef TPCSCGL_H
#define TPCSCGL_H
struct tpcSCGL_st {
	float SC[8]; /* Scale factor relating luminosity scaler to SpaceCharge */
	float SCoffset[8]; /* Offset to define luminosity for SpaceCharge */
	float SCexponent[8]; /* Luminosity exponential factor for SpaceCharge */
	float SCscaler[8]; /* Luminosity detector scaler */
	float GL[24]; /* Scale factor relating SpaceCharge to GridLeak */
	float GLoffset[24]; /* Offset to define luminosity for GridLeak */
	float GLradius; /* Radius of GridLeak between inner/outer sectors */
	float GLwidth; /* Width of GridLeak between inner/outer sectors */
	int mode; /* Modes to simplify parameter controls */
	char comment[256]; 
};
#endif /* TPCSCGL_H */
