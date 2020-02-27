/* TpcEffectivedX.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
:Description: Effective height of pad row
 */
#ifndef TPCEFFECTIVEDX_H
#define TPCEFFECTIVEDX_H
#define TPCEFFECTIVEDX_SPEC \
"struct TpcEffectivedX { \
	float scaleInner; \
	float scaleOuter; \
};"
typedef struct TpcEffectivedX_st {
	float scaleInner; /* scale factor for inner dX */
	float scaleOuter; /*       -"-        outer dX  */
} TPCEFFECTIVEDX_ST;
#endif /* TPCEFFECTIVEDX_H */
