/* trigDetSums.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    trigDetSums.idl

  Table: trigDetSums

       description:

 */
#ifndef TRIGDETSUMS_H
#define TRIGDETSUMS_H
#define TRIGDETSUMS_SPEC \
"struct trigDetSums { \
	unsigned long runNumber; \
	unsigned long timeOffset; \
	double ctbWest; \
	double ctbEast; \
	double ctbTOFp; \
	double tofp; \
	double zdcWest; \
	double zdcEast; \
	double zdcX; \
	double mult; \
	double L0; \
	double bbcX; \
	double bbcXctbTOFp; \
	double bbcWest; \
	double bbcEast; \
	double bbcYellowBkg; \
	double bbcBlueBkg; \
	double pvpdWest; \
	double pvpdEast; \
};"
typedef struct trigDetSums_st {
	unsigned int runNumber; /*       run number  */
	unsigned int timeOffset; /*       run begin time  */
	double ctbWest; /*   ctb West  */
	double ctbEast; /*   ctb East  */
	double ctbTOFp; /*   ctbOr + TOFp rate  */
	double tofp; /*   TOFp rate  */
	double zdcWest; /*    zdc west rate  */
	double zdcEast; /*    zdc east rate  */
	double zdcX; /*   zdc and rate  */
	double mult; /*   mult rate  */
	double L0; /*   L0 Rate  */
	double bbcX; /*   BBC and Rate  */
	double bbcXctbTOFp; /*   BBCAnd + ctbTOFp rate  */
	double bbcWest; /*   --BBC West--  */
	double bbcEast; /*   --BBC East--  */
	double bbcYellowBkg; /*   --(BBC Eastdelayed) and (BBC West)--  */
	double bbcBlueBkg; /*   --(BBC Westdelayed) and (BBC East)--  */
	double pvpdWest; /*   --PVPD East--  */
	double pvpdEast; /*   --PVPD West--  */
} TRIGDETSUMS_ST;
#endif /* TRIGDETSUMS_H */
