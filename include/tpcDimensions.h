/* tpcDimensions.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcDimensions.idl

  Table: tpcDimensions

       description:

 */
#ifndef TPCDIMENSIONS_H
#define TPCDIMENSIONS_H
struct tpcDimensions_st {
  int numberOfSectors; /*   */
  double tpcInnerRadius; /*   */
  double tpcOuterRadius; /*   */
  double tpcTotalLength; /*   */
  double wheelInnerRadius; /*   */
  double wheelOuterRadius; /*   */
  double wheelThickness; /*   */
  double senseGasOuterRadius; /*   */
  double tpeaThickness; /*   */
  double cathodeInnerRadius; /*   */
  double cathodeOuterRadius; /*   */
  double cathodeThickness; /*   */
  double outerCuThickness; /*   */
  double outerKaptonThickness; /*   */
  double outerNomexThickness; /*   */
  double outerGlueThickness; /*   */
  double outerInsGasThickness; /*   */
  double outerAlThickness; /*   */
  double outerAlHoneycombThickness; /*   */
  double innerGlueThickness; /*   */
  double innerNomexThickness; /*   */
  double innerKaptonThickness; /*   */
  double innerAlThickness; /*   */
  double innerGapWidI; /*  inner width of air in support wheel  */
  double innerGapWidO; /*  outer width of air in support wheel  */
  double innerGapHeit; /*  height (dr) of air in support wheel  */
  double innerGapRad; /*  air in support wheel - center radius  */
  double innerInWidth; /*   sector width at inner radius  */
  double innerOutWidth; /*   sector width at outer radius  */
  double innerHeight; /*   sector radial height  */
  double innerPPDepth; /*  padplane thickness (Al+pcb)  */
  double innerAlDepth; /*  depth of openings in Al structure  */
  double innerMWCDepth; /*  full thickness MWC sensitive region  */
  double innerBoundary; /*  Al frame boundary width  */
  double innerRCenter; /*  sector center radius (precision hole)  */
  double innerMWCInn; /*  MWC sensitive region inner size  */
  double innerMWCOut; /*  MWC sensitive region outer size  */
  double innerMVCHei; /*  MVC sensitive region radial  */
  int innerAirGaps; /*   number of air gaps in Al  */
  int innerExtraAl; /*   number of extra Al supportpieces  */
  double innerZGaps[5]; /*   opening positions  */
  double innerZGapsSize[5]; /*   opening size  */
  double innerXExtraAl[5]; /*   x positions  */
  double innerZExtraAl[5]; /*   z positions  */
  double innerDXExtraAl[5]; /*   x thickness  */
  double innerDZExtraAl[5]; /*   z thickness  */
  double outerGapWidI; /*  inner width of air in support wheel  */
  double outerGapWidO; /*  outer width of air in support wheel  */
  double outerGapHeit; /*  height (dr) of air in support wheel  */
  double outerGapRad; /*  air in support wheel - center radius  */
  double outerInWidth; /*   sector width at inner radius  */
  double outerOutWidth; /*   sector width at outer radius  */
  double outerHeight; /*   sector radial height  */
  double outerPPDepth; /*  padplane thickness (Al+pcb)  */
  double outerAlDepth; /*  depth of openings in Al structure  */
  double outerMWCDepth; /*  full thickness MWC sensitive region  */
  double outerBoundary; /*  Al frame boundary width  */
  double outerRCenter; /*  sector center radius (precision hole)  */
  double outerMWCInn; /*  MWC sensitive region inner size  */
  double outerMWCOut; /*  MWC sensitive region outer size  */
  double outerMVCHei; /*  MVC sensitive region radial  */
  int outerAirGaps; /*   number of air gaps in Al  */
  int outerExtraAl; /*   number of extra Al supportpieces  */
  double outerZGaps[8]; /*   opening positions  */
  double outerZGapsSize[8]; /*   opening size  */
  double outerXExtraAl[5]; /*   x positions  */
  double outerZExtraAl[5]; /*   z positions  */
  double outerDXExtraAl[5]; /*   x thickness  */
  double outerDZExtraAl[5]; /*   z thickness  */
};
#endif /* TPCDIMENSIONS_H */
