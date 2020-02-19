/* tpcWirePlanes.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcWirePlanes.idl

  Table: tpcWirePlanes

       description:

 */
#ifndef TPCWIREPLANES_H
#define TPCWIREPLANES_H
struct tpcWirePlanes_st {
  double anodeWireRadius; /*   */
  double frischGridWireRadius; /*   */
  double gatingGridWireRadius; /*   */
  double anodeWirePitch; /*   */
  double frischGridWirePitch; /*   */
  double gatingGridWirePitch; /*   */
  double innerSectorAnodeWirePadSep; /*   AnodeWire-to-PadPlane distance  */
  double innerSectorFrischGridPadSep; /*   FrischGrid-to-PadPlane distance  */
  double innerSectorGatingGridPadSep; /*   GatingGrid-to-PadPlane distance  */
  double outerSectorAnodeWirePadSep; /*   AnodeWire-to-PadPlane distance  */
  double outerSectorFrischGridPadSep; /*   FrischGrid-to-PadPlane distance  */
  double outerSectorGatingGridPadSep; /*   GatingGrid-to-PadPlane distance  */
  int numInnerSectorAnodeWires; /*   */
  int numInnerSectorFrischGridWires; /*   */
  int numInnerSectorGatingGridWires; /*   */
  double firstInnerSectorAnodeWire; /*   */
  double firstInnerSectorFrischGridWire; /*   */
  double firstInnerSectorGatingGridWire; /*   */
  double lastInnerSectorAnodeWire; /*   */
  int numOuterSectorAnodeWires; /*   */
  int numOuterSectorFrischGridWires; /*   */
  int numOuterSectorGatingGridWires; /*   */
  double firstOuterSectorAnodeWire; /*   */
  double firstOuterSectorFrischGridWire; /*   */
  double firstOuterSectorGatingGridWire; /*   */
  double lastOuterSectorAnodeWire; /*   */
};
#endif /* TPCWIREPLANES_H */
