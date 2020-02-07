/* tpcEffectiveGeom.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcEffectiveGeom.idl

  Table: tpcEffectiveGeom

       description:

 */
#ifndef TPCEFFECTIVEGEOM_H
#define TPCEFFECTIVEGEOM_H
#define TPCEFFECTIVEGEOM_SPEC \
"struct tpcEffectiveGeom { \
	double drift_length_correction; \
	double z_inner_offset; \
	double z_outer_offset; \
	double z_inner_offset_West; \
	double z_outer_offset_West; \
};"
typedef struct tpcEffectiveGeom_st {
  double drift_length_correction; /*  cm: Diff between actual drift length and  */
  double z_inner_offset; /*  cm: Effective distance between  */
  double z_outer_offset; /*  cm: Effective distance between  */
  double z_inner_offset_West; /*  cm: Effective distance West with respect to East */
  double z_outer_offset_West; /*  cm: Effective distance West  -"-                 */
} TPCEFFECTIVEGEOM_ST;
#endif /* TPCEFFECTIVEGEOM_H */
