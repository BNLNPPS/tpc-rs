/* tpcElectronics.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcElectronics.idl

  Table: tpcElectronics

       description:

 */
#ifndef TPCELECTRONICS_H
#define TPCELECTRONICS_H
#define TPCELECTRONICS_SPEC \
"struct tpcElectronics { \
	long numberOfTimeBins; \
	double nominalGain; \
	double samplingFrequency; \
	double tZero; \
	double adcCharge; \
	double adcConversion; \
	double averagePedestal; \
	double shapingTime; \
	double tau; \
};"
typedef struct tpcElectronics_st {
  int numberOfTimeBins; /*   */
  double nominalGain; /*   mV/fC  */
  double samplingFrequency; /* MHz, not used,  overwritten by starClockOnl*/
  double tZero; /*   us (microseconds)  */
  double adcCharge; /*   fC/adc count  */
  double adcConversion; /*   mV/adc count  */
  double averagePedestal; /*   adc counts  */
  double shapingTime; /*   ns  */
  double tau; /*   ns  */
} TPCELECTRONICS_ST;
#endif /* TPCELECTRONICS_H */
