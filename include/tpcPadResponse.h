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
struct tpcPadResponse_st {
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
};
#endif /* TPCPADRESPONSE_H */
