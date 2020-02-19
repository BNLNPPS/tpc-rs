/* tpcHighVoltages.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcHighVoltages.idl

  Table: tpcHighVoltages

       description:   Cathode and gating grid voltage for ExB distortions  Cathode and gating grid voltage for ExB distortions


 */
#ifndef TPCHIGHVOLTAGES_H
#define TPCHIGHVOLTAGES_H
struct tpcHighVoltages_st {
	float cathode; /*   kVolts  */
	float gatedGridRef; /*   Volts - nominal TPC value but is set by 48 sub-sectors  */
	float gridLeakWallTip[24]; /*   Volts - iTPC GridLeak wall tip voltage for 24 sectors  */
	float gridLeakWallSide[24]; /*   above +100 means no wall  */
};
#endif /* TPCHIGHVOLTAGES_H */
