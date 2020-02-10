/* tpcCorrection.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
:Description: Drift Distance depended correction
:Synonyms::::
:Source:
:Update:
:Update frequncy:
:Reminder:
:Recall frequency:
:Size of Data:
:Pointer to data:  tpcCorrection.time:
 type = 0 polymonical fit,                                        use only [min,max]
 type = 1 TChebyshev poly in range [min,max] => [-1,1]           
 type = 2 shifted  TChebyshev poly in range [min,max] => [ 0,1]  
 type = 3 X => Log(1 - |x|)                                       use only [min,max]
 type = 4 X => Log(1 - |x|)*sign(x)                                -"-
 type = 5 X => Log(x), for x <= 1 => 0
 type = 10 (== log(1. + OffSet/x) + poly(x,npar))                  -"-
 type = 11 (== log(1. + OffSet/x) + poly(x,npar) for log(ADC) and |Z|
 type = 200 cut on range [min,max]
 type = 300 don't correct out of range [min,max]
 type = 1000      ; gaus(0)+pol0(3);
 type = 1000 + 100; gaus(0)+pol1(3)
 type = 1000 + 200; gaus(0)+pol2(3)
 type = 1000 + 300; gaus(0)+pol3(3)
 type = 2000      ; expo(0)+pol0(2);
 type = 2000 + 100; expo(0)+pol1(2)
 type = 2000 + 200; expo(0)+pol2(2)
 type = 2000 + 300; expo(0)+pol3(2)
 row index

COMMENTS TRUNCATED */
#ifndef TPCCORRECTION_H
#define TPCCORRECTION_H
#define TPCCORRECTION_SPEC \
"struct tpcCorrection { \
	long type; \
	long idx; \
	long nrows; \
	long npar; \
	double OffSet; \
	double min; \
	double max; \
	double a[10]; \
};"
typedef struct tpcCorrection_st {
	int type; 
	int idx; 
	int nrows; 
	int npar; 
	double OffSet; 
	double min; 
	double max; 
	double a[10]; 
} TPCCORRECTION_ST;
#endif /* TPCCORRECTION_H */
