/* tpcPadResponse.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcPadResponse.idl

  Table: tpcPadResponse

       description:   trs parameters for pad response functions  trs parameters for pad response functions


 */
#ifndef TPCPADRESPONSE_H
#define TPCPADRESPONSE_H
#define TPCPADRESPONSE_SPEC \
"struct tpcPadResponse { \
	float innerGasGainFluctuation; \
	float outerGasGainFluctuation; \
	float innerPadResponseSigma; \
	float outerPadResponseSigma; \
	float innerWirePadCoupling; \
	float outerWirePadCoupling; \
	float innerRowNormalization; \
	float outerRowNormalization; \
	float BoundaryOfStepFunctions[6]; \
	float innerChargeFractionConstants[6]; \
	float outerChargeFractionConstants[6]; \
	float errorFunctionRange; \
	long errorFunctionEntry; \
	float longitudinalDiffusionConstant; \
	float transverseDiffusionConstant; \
	float InnerOuterFactor; \
};"
typedef struct tpcPadResponse_st {
  float innerGasGainFluctuation; /*   unitless  */
  float outerGasGainFluctuation; /*   unitless  */
  float innerPadResponseSigma; /*   cm  */
  float outerPadResponseSigma; /*   cm  */
  float innerWirePadCoupling; /*   cm  */
  float outerWirePadCoupling; /*   cm  */
  float innerRowNormalization; /*   unitless  */
  float outerRowNormalization; /*   unitless  */
  float BoundaryOfStepFunctions[6]; /*   cm  */
  float innerChargeFractionConstants[6]; /*   unitless  */
  float outerChargeFractionConstants[6]; /*   unitless  */
  float errorFunctionRange; /*   unitless  */
  int errorFunctionEntry; /*   unitless  */
  float longitudinalDiffusionConstant; /*   cm/sqrt(cm)  */
  float transverseDiffusionConstant; /*   cm/sqrt(cm)  */
  float InnerOuterFactor; /*   dimensionless  */
} TPCPADRESPONSE_ST;
#endif /* TPCPADRESPONSE_H */
