/* tpcRDOMasks.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcRDOMasks.idl

  Table: tpcRDOMasks

       description:

 */
#ifndef TPCRDOMASKS_H
#define TPCRDOMASKS_H
#define TPCRDOMASKS_SPEC \
"struct tpcRDOMasks { \
	unsigned long runNumber; \
	unsigned long sector; \
	unsigned long mask; \
};"
typedef struct tpcRDOMasks_st {
	unsigned int runNumber; /*       run number  */
	unsigned int sector; /*   sector  */
	unsigned int mask; /*   enable mask  */
} TPCRDOMASKS_ST;
#endif /* TPCRDOMASKS_H */
