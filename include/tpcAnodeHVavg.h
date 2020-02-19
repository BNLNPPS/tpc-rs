/* tpcAnodeHVavg.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcAnodeHV.idl

  Table: tpcAnodeHVavg

       description: average voltages over stable periods of time,

 */
#ifndef TPCANODEHVAVG_H
#define TPCANODEHVAVG_H
struct tpcAnodeHVavg_st {
  unsigned short sector; /*  sector 1-24 */
  unsigned short socket; /*  MWC socket/card (ISOR=17,OSIR=18,OSOR=19)  */
  float voltage; /*  average voltage  */
  float rms; /*  rms for averaged voltage */
  int numentries; /*  number of entries used for average */
  int numoutliers; /*  number of encountered outliers */
};
#endif /* TPCANODEHVAVG_H */
