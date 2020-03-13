/* beamInfo.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    beamInfo.idl

  Table: beamInfo

       description:     beam state for run at start and end   beam state for run at start and end


 */
#ifndef BEAMINFO_H
#define BEAMINFO_H
#define BEAMINFO_SPEC \
"struct beamInfo { \
	unsigned long runNumber; \
	long entryTag; \
	char blueSpecies[32]; \
	unsigned long blueMassNumber; \
	float blueEnergy; \
	float blueIntensity; \
	float blueLifeTime; \
	float blueBunchIntensity; \
	char yellowSpecies[32]; \
	unsigned long yellowMassNumber; \
	float yellowEnergy; \
	float yellowIntensity; \
	float yellowLifeTime; \
	float yellowBunchIntensity; \
	float blueFillNumber; \
	float yellowFillNumber; \
};"
typedef struct beamInfo_st {
	unsigned int runNumber; /*   */
	int entryTag; /*     0=startrun, 1=endrun, 2=runave, 3=std  */
	char blueSpecies[32]; /*     species  */
	unsigned int blueMassNumber; /*   */
	float blueEnergy; /*     energy  */
	float blueIntensity; /*     Ions  */
	float blueLifeTime; /*     Ions per minute  */
	float blueBunchIntensity; /*     bunch intensity  */
	char yellowSpecies[32]; /*     species  */
	unsigned int yellowMassNumber; /*   */
	float yellowEnergy; /*     energy  */
	float yellowIntensity; /*     Ions  */
	float yellowLifeTime; /*     Ions per minute  */
	float yellowBunchIntensity; /*     bunch intensity  */
	float blueFillNumber; /*   */
	float yellowFillNumber; /*   */
} BEAMINFO_ST;
#endif /* BEAMINFO_H */