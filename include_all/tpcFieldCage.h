/* tpcFieldCage.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
:Description: 
  :Synonyms::::
  :Source:
  :Update:
  :Update frequncy:
  :Reminder:
  :Recall frequency:
  :Size of Data:
  :Pointer to data:  

 */
#ifndef TPCFIELDCAGE_H
#define TPCFIELDCAGE_H
#define TPCFIELDCAGE_SPEC \
"struct tpcFieldCage { \
	float innerFieldCageShift; \
	float eastClockError; \
	float westClockError; \
};"
typedef struct tpcFieldCage_st {
	float innerFieldCageShift; /* cm : z shift of inner field cage w.r.t outer field cage */
	float eastClockError; /* radians :  Phi rotation of East end of TPC in radians */
	float westClockError; /* radians :  Phi rotation of West end of TPC in radians */
} TPCFIELDCAGE_ST;
#endif /* TPCFIELDCAGE_H */
