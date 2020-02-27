/* TpcSecRowCor.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    TpcSecRowCor.idl

  Table: TpcSecRowCor
  
       description: Tpc gain correction for sector  row

 */
#ifndef TPCSECROWCOR_H
#define TPCSECROWCOR_H
#define TPCSECROWCOR_SPEC \
"struct TpcSecRowCor { \
	float GainScale[100]; \
	float GainRms[100]; \
};"
typedef struct TpcSecRowCor_st {
	float GainScale[100]; /*  Gains for sector & row */
	float GainRms[100]; /*  RMS  - " -  */
} TPCSECROWCOR_ST;
#endif /* TPCSECROWCOR_H */
