/* Survey.h */
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
:Pointer to data:  Survey.time.C:
   Translation from local to Master
   m = R*l + t
           ( r00 r01 r02 ) (xl)   ( t0 )    ( r00 r01 r02 ) (xl)   ( t0 )    (     1 -gamma   beta ) (xl)   ( t0 )
       R = ( r10 r11 r12 ) (yl) + ( t1 ) =  ( r10 r11 r12 ) (zl) + ( t1 ) ~  ( gamma      1 -alpha ) (zl) + ( t1 ) 
           ( r20 r21 r22 ) (zl)   ( t2 )    ( r20 r21 r22 ) (yl)   ( t2 )    ( -beta  alpha      1 ) (yl)   ( t2 )
 SVT
 Id = 0                                for SvtOnGlobal
 Id = [0, 1]                           for ShellOnGlobal, 
 0 is the x (South) Shell, 1 is the -x (North) Shell"
 Id = 1000*barrel + ladder             for LadderOnSurvey
 Id = 1000*barrel + ladder             for LadderOnShell
 Id = 1000*barrel + 100*wafer + ladder for WaferOnLadder

 SSD
 Id = 0                                for SsdOnGlobal
 Id = sector [1-4]                     SsdSectorsOnGlobal
 Id = 100*sector + ladder              SsdLaddersOnSectors
 Id = 7000 + 100*wafer + ladder        SsdWafersOnLadders
 */
#ifndef SURVEY_H
#define SURVEY_H
struct Survey_st {
	int Id; 
	double r00; 
	double r01; /* -gamma */
	double r02; /*  beta  */
	double r10; /*  gamma */
	double r11; 
	double r12; /* -alpha */
	double r20; /* -beta  */
	double r21; /*  alpha */
	double r22; 
	double t0; 
	double t1; 
	double t2; 
	double sigmaRotX; 
	double sigmaRotY; 
	double sigmaRotZ; 
	double sigmaTrX; 
	double sigmaTrY; 
	double sigmaTrZ; 
	char comment[32]; 
};
#endif /* SURVEY_H */
