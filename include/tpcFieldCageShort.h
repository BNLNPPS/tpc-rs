/* tpcFieldCageShort.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcFieldCageShort.idl

  Table: tpcFieldCageShort

       description:     provide info on shorted rings in the field cages

 */
#ifndef TPCFIELDCAGESHORT_H
#define TPCFIELDCAGESHORT_H
struct tpcFieldCageShort_st {
	float side; /* 0 = east, 1 = west */
	float cage; /* 0 = inner, 1 = outer */
	float ring; /* ring location of the short (e.g. 169.5) */
	float resistor; /* MOhm value of added external resistor to resistor chain */
	float MissingResistance; /* missing resistance */
};
#endif /* TPCFIELDCAGESHORT_H */
