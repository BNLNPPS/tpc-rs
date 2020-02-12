/* tpcRDOMap.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcRDOMap.idl

  Table: tpcRDOMap

       description: Map tpc row, pad range [padMin, padMax] to rdo number

 row index
 total no. of real rows in the table; For Db interface (where nrows = 50)
 */
#ifndef TPCRDOMAP_H
#define TPCRDOMAP_H
#define TPCRDOMAP_SPEC \
"struct tpcRDOMap { \
	long idx; \
	long nrows; \
	octet row; \
	octet padMin; \
	octet padMax; \
	octet rdo; \
};"
typedef struct tpcRDOMap_st {
	int idx; 
	int nrows; 
	unsigned char row; 
	unsigned char padMin; 
	unsigned char padMax; 
	unsigned char rdo; 
} TPCRDOMAP_ST;
#endif /* TPCRDOMAP_H */
