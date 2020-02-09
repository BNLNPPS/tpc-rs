/* tpcAnodeHV.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcAnodeHV.idl

  Table: tpcAnodeHV

       description:

 */
#ifndef TPCANODEHV_H
#define TPCANODEHV_H
#define TPCANODEHV_SPEC \
"struct tpcAnodeHV { \
	unsigned short sector; \
	unsigned short socket; \
	float voltage; \
};"
typedef struct tpcAnodeHV_st {
	unsigned short sector; /*  sector 1-24 */
	unsigned short socket; /*  MWC socket/card (ISOR=17,OSIR=18,OSOR=19)  */
	float voltage; /*   HV setting  */
} TPCANODEHV_ST;
#endif /* TPCANODEHV_H */
