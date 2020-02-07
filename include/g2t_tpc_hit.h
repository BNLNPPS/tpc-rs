/* g2t_tpc_hit.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
#ifndef G2T_TPC_HIT_H
#define G2T_TPC_HIT_H
#define G2T_TPC_HIT_SPEC \
"struct g2t_tpc_hit { \
	long id; \
	long next_tr_hit_p; \
	long track_p; \
	long volume_id; \
	float de; \
	float ds; \
	float p[3]; \
	float tof; \
	float x[3]; \
	float lgam; \
	float length; \
	float adc; \
	float pad; \
	float timebucket; \
	long np; \
};"
typedef struct g2t_tpc_hit_st {
  int id; /* primary key */
  int next_tr_hit_p; /* Id of next hit on same track */
  int track_p; /* Id of parent track */
  int volume_id; /* STAR volume identification */
  float de; /* energy deposition at hit */
  float ds; /* path length within padrow */
  float p[3]; /* local momentum */
  float tof; /* time of flight */
  float x[3]; /* coordinate (Cartesian) */
  float lgam; /* ALOG10(GEKin/AMass) */
  float length; /* track length up to this hit */
  float adc; /* signal in ADC  after digitization */
  float pad; /* hit pad position used in digitization */
  float timebucket; /* hit time position -"- */
  int np; /* no. of primary electrons */
} G2T_TPC_HIT_ST;
#endif /* G2T_TPC_HIT_H */
