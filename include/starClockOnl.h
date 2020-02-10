/* starClockOnl.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    starClockOnl.idl

  Table: starClockOnl

       description:

 */
#ifndef STARCLOCKONL_H
#define STARCLOCKONL_H
#define STARCLOCKONL_SPEC \
"struct starClockOnl { \
	unsigned long runNumber; \
	unsigned long time; \
	double frequency; \
};"
typedef struct starClockOnl_st {
	unsigned int runNumber; /*   run number  */
	unsigned int time; /*   unix time of entry  */
	double frequency; /*   frequency in Hz  */
} STARCLOCKONL_ST;
#endif /* STARCLOCKONL_H */
