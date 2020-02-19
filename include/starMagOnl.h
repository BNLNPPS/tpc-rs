/* starMagOnl.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    starMagOnl.idl

  Table: starMagOnl

       description:

 */
#ifndef STARMAGONL_H
#define STARMAGONL_H
struct starMagOnl_st {
	unsigned int runNumber; /*   run number  */
	unsigned int time; /*   unix time of entry  */
	double current; /*   magnet current (- means B polarity)  */
};
#endif /* STARMAGONL_H */
