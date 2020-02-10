/* spaceChargeCor.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    spaceChargeCor.idl

  Table: spaceChargeCor

       description: Table for Space Charge Corrections

 */
#ifndef SPACECHARGECOR_H
#define SPACECHARGECOR_H
#define SPACECHARGECOR_SPEC \
"struct spaceChargeCor { \
	double fullFieldB; \
	double halfFieldB; \
	double zeroField; \
	double halfFieldA; \
	double fullFieldA; \
	double satRate; \
	float factor; \
	float detector; \
	float offset; \
	float ewratio; \
};"
typedef struct spaceChargeCor_st {
	double fullFieldB; /* Negative Full Field Correction  */
	double halfFieldB; /* Negative Half Field Correction  */
	double zeroField; /*  Zero Field " "  */
	double halfFieldA; /*  Postive Half " " */
	double fullFieldA; /*  Postive Full " " */
	double satRate; /* Saturation Rate Hz  */
	float factor; /*  Multiplicative Factor */
	float detector; /* 0=VPDx, 1=BBCx, 2=ZDCx, 3=ZDCe+w, 4=BBCe+w, ... */
	float offset; /* Offset at zero luminosity */
	float ewratio; /* Ratio of charge east/west */
} SPACECHARGECOR_ST;
#endif /* SPACECHARGECOR_H */
