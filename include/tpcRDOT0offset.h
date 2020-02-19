/* tpcRDOT0offset.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcRDOT0offset.idl

  Table: tpcRDOT0offset

       description:

 */
#ifndef TPCRDOT0OFFSET_H
#define TPCRDOT0OFFSET_H
struct tpcRDOT0offset_st {
	unsigned char isShifted[24]; /* flag if there is any RDO off set to tsector */
	float t0[24][10]; /* RDO t0 offset per sector: [0-23], rdo [0-5] Tpx, [6-9] iTpc  (time bins) */
};
#endif /* TPCRDOT0OFFSET_H */
