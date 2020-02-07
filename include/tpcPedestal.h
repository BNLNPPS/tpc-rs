/* tpcPedestal.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcPedestal.idl

  Table: tpcPedestal

       description:

 */
#ifndef TPCPEDESTAL_H
#define TPCPEDESTAL_H
#define TPCPEDESTAL_SPEC \
"struct tpcPedestal { \
	float Pedestal[100][182]; \
	float Rms[100][182]; \
};"
typedef struct tpcPedestal_st {
  float Pedestal[100][182]; /*  Pedestals per */
  float Rms[100][182]; /*  Rms per */
} TPCPEDESTAL_ST;
#endif /* TPCPEDESTAL_H */
