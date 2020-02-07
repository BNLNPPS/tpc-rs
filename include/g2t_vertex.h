/* g2t_vertex.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
#ifndef G2T_VERTEX_H
#define G2T_VERTEX_H
#define G2T_VERTEX_SPEC \
"struct g2t_vertex { \
	char ge_volume[4]; \
	long daughter_p; \
	long eg_label; \
	long eg_proc; \
	long event_p; \
	long ge_medium; \
	long ge_proc; \
	long id; \
	long is_itrmd; \
	long n_daughter; \
	long n_parent; \
	long next_itrmd_p; \
	long next_prim_v_p; \
	long parent_p; \
	float eg_tof; \
	float eg_x[3]; \
	float ge_tof; \
	float ge_x[3]; \
};"
typedef struct g2t_vertex_st {
  char ge_volume[4]; /* GEANT volume name */
  int daughter_p; /* Id of first daughter in linked list */
  int eg_label; /* generator label (or 0 if GEANT vertex) */
  int eg_proc; /* generator production mechanism (if any) */
  int event_p; /* pointer to event */
  int ge_medium; /* GEANT Medium */
  int ge_proc; /* >0 GEANT, =0 event generator */
  int id; /* primary key */
  int is_itrmd; /* flags intermediate vertex */
  int n_daughter; /* Number of daughter tracks */
  int n_parent; /* number of parent tracks */
  int next_itrmd_p; /* Id of next intermedate vertex */
  int next_prim_v_p; /* Id of next primary vertex */
  int parent_p; /* Id of first parent track */
  float eg_tof; /* generator vertex production time */
  float eg_x[3]; /* generator vertex coordinate (Cartesian) */
  float ge_tof; /* GEANT vertex production time */
  float ge_x[3]; /* GEANT vertex coordinate (Cartesian) */
} G2T_VERTEX_ST;
#endif /* G2T_VERTEX_H */
