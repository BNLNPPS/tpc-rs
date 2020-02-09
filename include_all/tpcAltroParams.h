/* tpcAltroParams.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcAltroParams.idl
 
   Table: tpcAltroParams
 
        description: Parameters of Altro cheap configuration
        N  no. of Altro parameters, 
        N  0  > filter is switched off
        N  1 > old TPC electronics
        Altro::ConfigZerosuppression(int Threshold, int MinSamplesaboveThreshold, int Presamples, int Postsamples)
        Altro::ConfigTailCancellationFilter(int K1, int K2, int K3, int L1, int L2, int L3)
 
  Threshold
  MinSamplesaboveThreshold
  K1 coefficient of the TCF
  K2 coefficient of the TCF
  K3 coefficient of the TCF
  L1 coefficient of the TCF
  L2 coefficient of the TCF
  L3 coefficient of the TCF
 */
#ifndef TPCALTROPARAMS_H
#define TPCALTROPARAMS_H
#define TPCALTROPARAMS_SPEC \
"struct tpcAltroParams { \
	long N; \
	long Altro_thr; \
	long Altro_seq; \
	long Altro_K1; \
	long Altro_K2; \
	long Altro_K3; \
	long Altro_L1; \
	long Altro_L2; \
	long Altro_L3; \
};"
typedef struct tpcAltroParams_st {
	int N; /*  = no. of Altro parameters */
	int Altro_thr; 
	int Altro_seq; 
	int Altro_K1; 
	int Altro_K2; 
	int Altro_K3; 
	int Altro_L1; 
	int Altro_L2; 
	int Altro_L3; 
} TPCALTROPARAMS_ST;
#endif /* TPCALTROPARAMS_H */
