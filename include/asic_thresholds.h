/* asic_thresholds.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
#ifndef ASIC_THRESHOLDS_H
#define ASIC_THRESHOLDS_H
#define ASIC_THRESHOLDS_SPEC \
"struct asic_thresholds { \
	long thresh_lo; \
	long thresh_hi; \
	long n_seq_lo; \
	long n_seq_hi; \
};"
typedef struct asic_thresholds_st {
  int thresh_lo;
  int thresh_hi;
  int n_seq_lo;
  int n_seq_hi;
} ASIC_THRESHOLDS_ST;
#endif /* ASIC_THRESHOLDS_H */
