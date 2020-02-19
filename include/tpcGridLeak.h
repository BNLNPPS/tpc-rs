/* tpcGridLeak.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcGridLeak.idl

  Table: tpcGridLeak

       description: Table for TPC GridLeaks

 */
#ifndef TPCGRIDLEAK_H
#define TPCGRIDLEAK_H
struct tpcGridLeak_st {
	double InnerGLRadius; /* Radius of GL around inner sectors           */
	double MiddlGLRadius; /* Radius of GL between inner/outer sectors    */
	double OuterGLRadius; /* Radius of GL around outer sectors           */
	double InnerGLWidth; /* Width of GL around inner sectors            */
	double MiddlGLWidth; /* Width of GL between inner/outer sectors     */
	double OuterGLWidth; /* Width of GL around outer sectors            */
	double InnerGLStrength; /* Strength of GL around inner sectors         */
	double MiddlGLStrength; /* Strength of GL between inner/outer sectors  */
	double OuterGLStrength; /* Strength of GL around outer sectors         */
};
#endif /* TPCGRIDLEAK_H */
