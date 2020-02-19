/* MDFCorrection.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
:Description: TMultiDimFit
:Synonyms::::
:Source:
:Update:
:Update frequncy:
:Reminder:
:Recall frequency:
:Size of Data:
:Pointer to data:  MDFCorrection.time:
 row index
 total no. of real rows in the table; For Db interface (where nrows = 50)
 type = 0 kMonomials, type = 1 kChebyshev, type = 2 kLegendre
 == 2 for now.
  p_ij = Power[i * NVariables + j];
 
 */
#ifndef MDFCORRECTION_H
#define MDFCORRECTION_H
struct MDFCorrection_st {
	unsigned char idx; 
	unsigned char nrows; 
	unsigned char PolyType; 
	unsigned char NVariables; 
	unsigned char NCoefficients; 
	unsigned char Power[100]; 
	double DMean; 
	double XMin[2]; 
	double XMax[2]; 
	double Coefficients[50]; 
	double CoefficientsRMS[50]; 
};
#endif /* MDFCORRECTION_H */
