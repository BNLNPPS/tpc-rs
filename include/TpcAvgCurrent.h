/* TpcAvgCurrent.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
#ifndef TPCAVGCURRENT_H
#define TPCAVGCURRENT_H
struct TpcAvgCurrent_st {
  int run; /* run no. used for averaging  */
  int start_time; /* begin unix time of averaging interval */
  int stop_time; /* end   unix time of averaging interval */
  float AvCurrent[192]; /* average current per sector(24) and channel(8) [muA]*/
  float AcCharge[192]; /* accumulated charge per sector(24) and channel(8) [C]*/
};
#endif /* TPCAVGCURRENT_H */
