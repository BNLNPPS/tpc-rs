/* tpcDriftVelocity.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcDriftVelocity.idl

  Table: tpcDriftVelocity

       description:

 */
#ifndef TPCDRIFTVELOCITY_H
#define TPCDRIFTVELOCITY_H
#define TPCDRIFTVELOCITY_SPEC \
"struct tpcDriftVelocity { \
	float laserDriftVelocityEast; \
	float laserDriftVelocityWest; \
	float cathodeDriftVelocityEast; \
	float cathodeDriftVelocityWest; \
};"
typedef struct tpcDriftVelocity_st {
	float laserDriftVelocityEast; /*   cm/us : from laser beam analysis  */
	float laserDriftVelocityWest; /*   cm/us : from laser beam analysis  */
	float cathodeDriftVelocityEast; /*   cm/us : from cathode emission  */
	float cathodeDriftVelocityWest; /*   cm/us : from cathode emission  */
} TPCDRIFTVELOCITY_ST;
#endif /* TPCDRIFTVELOCITY_H */
