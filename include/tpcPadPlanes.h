/* tpcPadPlanes.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcPadPlanes.idl

  Table: tpcPadPlanes

       description:

 */
#ifndef TPCPADPLANES_H
#define TPCPADPLANES_H
struct tpcPadPlanes_st {
  int padRows; /*   */
  int innerPadRows; /*   */
  int innerPadRows48; /*   */
  int innerPadRows52; /*   */
  int outerPadRows; /*   */
  int superInnerPadRows; /*   */
  int superOuterPadRows; /*   */
  double innerSectorPadWidth; /*   */
  double innerSectorPadLength; /*   */
  double innerSectorPadPitch; /*   */
  double innerSectorRowPitch1; /*   */
  double innerSectorRowPitch2; /*   */
  double firstPadRow; /*   */
  double firstOuterSectorPadRow; /*   */
  double lastOuterSectorPadRow; /*   */
  double firstRowWidth; /*   */
  double lastRowWidth; /*   */
  double outerSectorPadWidth; /*   */
  double outerSectorPadLength; /*   */
  double outerSectorPadPitch; /*   */
  double outerSectorRowPitch; /*   */
  double outerSectorLength; /*   */
  double ioSectorSeparation; /*   */
  double innerSectorEdge; /*   */
  double outerSectorEdge; /*   */
  double innerSectorPadPlaneZ; /*   */
  double outerSectorPadPlaneZ; /*   */
  int innerPadsPerRow[13]; /*   */
  int outerPadsPerRow[32]; /*   */
  double innerRowRadii[13]; /*   */
  double outerRowRadii[32]; /*   */
};
#endif /* TPCPADPLANES_H */
