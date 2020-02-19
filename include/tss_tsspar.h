/* tss_tsspar.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
   tss_tsspar_st.h
  Table containg parameters for slow simulator running.

 */
#ifndef TSS_TSSPAR_H
#define TSS_TSSPAR_H
struct tss_tsspar_st {
  char fileout[80]; /* output file for pixel data (none->table) */
  int dynam; /* adc dynamic range (adc counts; usu 1023) */
  int format; /* pixel data format */
  int max_itime; /* upper bound of time bucket */
  int max_pads; /* upper bound of pads */
  int max_row; /* upper bound of row (<=45) */
  int max_sect; /* upper bound of sector pair (<=12) */
  int min_itime; /* low bound of time bucket */
  int min_pads; /* lower bound of pads */
  int min_row; /* lower bound of row (>=1) */
  int min_sect; /* lower bound of sector pair (>=1) */
  int mode; /* mode of TPC simulation for diff. tests */
  int nele_laser; /* number of electrons from laser point*/
  int ngain; /* number of gain sampling (1-10) */
  int nseg; /* number of sub-segment within a G-vol */
  int ntime; /* number of time buckets (dimensionless) */
  int printout; /* control the level of printout (0=no) */
  int tpc_half; /* half (1) or full (0) TPC volume */
  int reset; /* re-do setup: 0=no, anything else=yes */
  float ave_ion_pot; /* Average Ion. Potential of a gas(Ar=26eV) */
  float bfield; /* magnetic field strength (Tesla) */
  float c_test; /* test capacitance value (pF) */
  float diff_long; /* long diff const of gas (cm/sqrt(cm)) */
  float diff_trans; /* trans diff const of gas (cm/sqrt(cm)) */
  float gain_in; /* gas gain:inner sector (dimensionless) */
  float gain_out; /* gas gain:outer sector (dimensionless) */
  float prf_in; /* pad resp.func:inner sector (cm) */
  float prf_out; /* pad resp.func:outer sector (cm) */
  float sca_rms; /* SCA noise (random, not filtered) */
  float scale; /* number of electrons per ADC count (600) */
  float step_size; /* step size for subpadrow tracking */
  float tau; /* shaper resp. time const. (usec) */
  float threshold; /* adc threshold (adc counts but a real num */
  float time_offset; /* Wayne Betts time offsets */
  float v_test; /* voltage of test pulse (V) */
  float white_rms; /* rms of white noise (shaper filtered) */
  float wire_coupling_in; /* wire-to-pad coupling inner sector */
  float wire_coupling_out; /* wire-to-pad coupling outer sector */
  float x_laser; /* local x of laser point[cm] along row */
  float y_laser; /* local y of laser point[cm] across row */
  float z_laser; /* z drift length of pointlaser source(cm) */
};
#endif /* TSS_TSSPAR_H */
