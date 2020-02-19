/* tpcPadGainT0.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcPadGainT0.idl

  Table: tpcPadGainT0

       description:

 */
#ifndef TPCPADGAINT0_H
#define TPCPADGAINT0_H
struct tpcPadGainT0_st {
  int run; /* pulser run number used */
  float Gain[24][45][182]; /* Gains per pad*/
  float T0[24][45][182]; /* T9 per pad*/
};
#endif /* TPCPADGAINT0_H */
