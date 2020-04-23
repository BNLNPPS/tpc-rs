#ifndef St_trigDetSumsC_h
#define St_trigDetSumsC_h

#include "tpcrs/config_structs.h"
#include <cmath>
#include "trigDetSums.h"
#include "StDetectorDbMaker/StDetectorDbClock.h"
#include "StDetectorDbMaker/St_richvoltagesC.h"
struct St_trigDetSumsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_trigDetSumsC, trigDetSums_st>
{
  unsigned int 	runNumber(int i = 0) 	        {return Struct(i)->runNumber;}
  unsigned int 	timeOffset(int i = 0) 	{return Struct(i)->timeOffset;}
  double 	ctbWest(int i = 0) 	        {return Struct(i)->ctbWest;}
  double 	ctbEast(int i = 0) 	        {return Struct(i)->ctbEast;}
  double 	ctbTOFp(int i = 0) 	        {return Struct(i)->ctbTOFp;}
  double 	tofp(int i = 0) 	        {return Struct(i)->tofp;}
  double 	zdcWest(int i = 0) 	        {return Struct(i)->zdcWest;}
  double 	zdcEast(int i = 0) 	        {return Struct(i)->zdcEast;}
  double 	zdcX(int i = 0) 	        {return Struct(i)->zdcX;}
  double 	mult(int i = 0) 	        {return Struct(i)->mult;}
  double 	L0(int i = 0) 	        {return Struct(i)->L0;}
  double 	bbcX(int i = 0) 	        {return Struct(i)->bbcX;}
  double 	bbcXctbTOFp(int i = 0) 	{return Struct(i)->bbcXctbTOFp;}
  double 	bbcWest(int i = 0) 	        {return Struct(i)->bbcWest;}
  double 	bbcEast(int i = 0) 	        {return Struct(i)->bbcEast;}
  double 	bbcYellowBkg(int i = 0) 	{return Struct(i)->bbcYellowBkg;}
  double 	bbcBlueBkg(int i = 0) 	{return Struct(i)->bbcBlueBkg;}
  double 	pvpdWest(int i = 0) 	        {return Struct(i)->pvpdWest;}
  double 	pvpdEast(int i = 0) 	        {return Struct(i)->pvpdEast;}
  double 	zdcCoin(int i = 0)            {return Nc(zdcX(i), zdcEast(i), zdcWest(i));}
  double 	bbcCoin(int i = 0)            {return Nc(bbcX(i), bbcEast(i), bbcWest(i));}
  void		validityMargin(double margin = 0) {fMargin = margin;}
  double getCTBWest() {return ctbWest();}
  double getCTBEast() {return ctbEast();}
  double getCTBOrTOFp() {return ctbTOFp();}
  double getTOFp() {return tofp();}
  double getZDCWest() {return zdcWest();}
  double getZDCEast() {return zdcEast();}
  double getZDCX() {return zdcX();}
  double getZDCCoin() {return zdcCoin();}
  double getMult() {return mult();}
  double getL0() {return L0();}
  double getBBCX() {return bbcX();}
  double getBBCCoin() {return bbcCoin();}
  double getBBCXCTB() {return bbcXctbTOFp();}
  double getBBCWest() {return bbcWest();}
  double getBBCEast() {return bbcEast();}
  double getBBCYellowBkg() {return bbcYellowBkg();}
  double getBBCBlueBkg() {return bbcBlueBkg();}
  double getPVPDWest() {return pvpdWest();}
  double getPVPDEast() {return pvpdEast();}
  unsigned int   getRichHVStatus() {return St_richvoltagesC::instance()->status();}
  void     setValidityMargin(double margin = 0) {validityMargin(margin);}

  // The following code attempts to correct coincidence rates for accidentals and multiples
  // See STAR Note 528
  static double Nc(double New, double Ne, double Nw, int n_bunches = 111)
  {
    // 111 is a guess using the maximum seen filled bunches in RHIC so far
    // (not always the case, but we don't have access to this number)
    double Nbc = StDetectorDbClock::instance()->CurrentFrequency() * ((double) n_bunches) / 120.;
    return -Nbc * std::log(1. - ((New - (Ne * Nw / Nbc)) / (Nbc + New - Ne - Nw)));
  }
 private:
  double	fMargin;
};
#endif
