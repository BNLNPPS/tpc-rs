/* tpcOmegaTau.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcOmegaTau.idl

  Table: tpcOmegaTau

       description: Table for TPC OmegaTau

 */
#ifndef TPCOMEGATAU_H
#define TPCOMEGATAU_H
struct tpcOmegaTau_st {
	float tensorV1; /* tensor for OmegaTau           */
	float tensorV2; /* tensor for OmegaTau    */
	unsigned short distortionCorrectionsMode; /* modes for field distortion and non-uniformity corrections  */
};
#endif /* TPCOMEGATAU_H */
