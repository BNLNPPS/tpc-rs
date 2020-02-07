/* TpcAvgPowerSupply.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
#ifndef TPCAVGPOWERSUPPLY_H
#define TPCAVGPOWERSUPPLY_H
#define TPCAVGPOWERSUPPLY_SPEC \
"struct TpcAvgPowerSupply { \
	long run; \
	long start_time; \
	long stop_time; \
	float Current[192]; \
	float Charge[192]; \
	float Voltage[192]; \
};"
typedef struct TpcAvgPowerSupply_st {
  int run; /* run no. used for averaging  */
  int start_time; /* begin unix time of averaging interval */
  int stop_time; /* end   unix time of averaging interval */
  float Current[192]; /* average current per sector(24) and channel(8) [muA]*/
  float Charge[192]; /* accumulated charge per sector(24) and channel(8) [C]*/
  float Voltage[192]; /* average Voltage per sector(24) and channel(8) [V]*/
} TPCAVGPOWERSUPPLY_ST;
#endif /* TPCAVGPOWERSUPPLY_H */
