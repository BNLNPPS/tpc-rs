/* tpcSectorT0offset.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcSectorT0offset.idl

  Table: tpcSectorT0offset

       description:

 */
#ifndef TPCSECTORT0OFFSET_H
#define TPCSECTORT0OFFSET_H
#define TPCSECTORT0OFFSET_SPEC \
"struct tpcSectorT0offset { \
	float t0[48]; \
};"
typedef struct tpcSectorT0offset_st {
	float t0[48]; /* Sector t0 offset per sector: [0-23] Tpx, [24-47] iTpc  (time bins)*/
} TPCSECTORT0OFFSET_ST;
#endif /* TPCSECTORT0OFFSET_H */
