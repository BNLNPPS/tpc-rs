/* tpcSlowControlSim.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcSlowControlSim.idl

  Table: tpcSlowControlSim

       description:

 */
#ifndef TPCSLOWCONTROLSIM_H
#define TPCSLOWCONTROLSIM_H
#define TPCSLOWCONTROLSIM_SPEC \
"struct tpcSlowControlSim { \
	double driftVelocity; \
	double driftVoltage; \
	double innerSectorAnodeVoltage; \
	double innerSectorGatingGridV; \
	double outerSectorAnodeVoltage; \
	double outerSectorGatingGridV; \
	double innerSectorGasGain; \
	double innerSectorGasGainVzero; \
	double innerSectorGasGainb; \
	double outerSectorGasGain; \
	double outerSectorGasGainVzero; \
	double outerSectorGasGainb; \
	double hallPressure; \
	double hallTemperature; \
};"
typedef struct tpcSlowControlSim_st {
  double driftVelocity; /*   */
  double driftVoltage; /*   */
  double innerSectorAnodeVoltage; /*   */
  double innerSectorGatingGridV; /*   */
  double outerSectorAnodeVoltage; /*   */
  double outerSectorGatingGridV; /*   */
  double innerSectorGasGain; /*   */
  double innerSectorGasGainVzero; /*   */
  double innerSectorGasGainb; /*   */
  double outerSectorGasGain; /*   */
  double outerSectorGasGainVzero; /*   */
  double outerSectorGasGainb; /*   */
  double hallPressure; /*   */
  double hallTemperature; /*   */
} TPCSLOWCONTROLSIM_ST;
#endif /* TPCSLOWCONTROLSIM_H */
