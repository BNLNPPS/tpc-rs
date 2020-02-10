/* iTPCSurvey.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
:Description: Survey data
:Synonyms::::
:Source:
:Update:
:Update frequncy:
:Reminder:
:Recall frequency:
:Size of Data:
:Pointer to data:  iTPCSurvey.time.C:
:iTPC survey measurements J.Thomas, 04/20/2018
 */
#ifndef ITPCSURVEY_H
#define ITPCSURVEY_H
#define ITPCSURVEY_SPEC \
"struct iTPCSurvey { \
	long Id; \
	float Angle; \
	float dx; \
	float dy; \
	float ScaleX; \
	float ScaleY; \
	char comment[32]; \
};"
typedef struct iTPCSurvey_st {
	int Id; 
	float Angle; 
	float dx; 
	float dy; 
	float ScaleX; 
	float ScaleY; 
	char comment[32]; 
} ITPCSURVEY_ST;
#endif /* ITPCSURVEY_H */
