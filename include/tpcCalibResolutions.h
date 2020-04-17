/* tpcCalibResolutions.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcCalibResolutions.idl

  Table: tpcCalibResolutions

       description: Table for Resolutions of TPC Calibrations

 */
#ifndef TPCCALIBRESOLUTIONS_H
#define TPCCALIBRESOLUTIONS_H
struct tpcCalibResolutions_st {
	float SpaceCharge; /* SpaceCharge correction */
	float GridLeak; /* GridLeak correction */
	char comment[255]; /* comments */
};
#endif /* TPCCALIBRESOLUTIONS_H */
