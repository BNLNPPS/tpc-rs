#ifndef St_trigDetSumsC_h
#define St_trigDetSumsC_h

#include "TChair.h"
#include <cmath>
#include "tables/St_trigDetSums_Table.h"
#include "StDetectorDbMaker/StDetectorDbClock.h"
#include "StDetectorDbMaker/St_richvoltagesC.h"
#include "TMath.h"
class St_trigDetSumsC : public TChair
{
 public:
  St_trigDetSumsC(St_trigDetSums* table = 0) : TChair(table) {SafeDelete(fgInstance); fgInstance = this; fMargin = 0;}
  virtual ~St_trigDetSumsC() {fgInstance = 0;}
  static St_trigDetSumsC* 	instance();
  static St_trigDetSumsC*      GetInstance() {return fgInstance;}
  trigDetSums_st* 	Struct(Int_t i = 0) 	{return ((St_trigDetSums*) Table())->GetTable() + i;}
  UInt_t     	getNumRows()                	{return GetNRows();}
  UInt_t 	runNumber(Int_t i = 0) 	        {return Struct(i)->runNumber;}
  UInt_t 	timeOffset(Int_t i = 0) 	{return Struct(i)->timeOffset;}
  Double_t 	ctbWest(Int_t i = 0) 	        {return Struct(i)->ctbWest;}
  Double_t 	ctbEast(Int_t i = 0) 	        {return Struct(i)->ctbEast;}
  Double_t 	ctbTOFp(Int_t i = 0) 	        {return Struct(i)->ctbTOFp;}
  Double_t 	tofp(Int_t i = 0) 	        {return Struct(i)->tofp;}
  Double_t 	zdcWest(Int_t i = 0) 	        {return Struct(i)->zdcWest;}
  Double_t 	zdcEast(Int_t i = 0) 	        {return Struct(i)->zdcEast;}
  Double_t 	zdcX(Int_t i = 0) 	        {return Struct(i)->zdcX;}
  Double_t 	mult(Int_t i = 0) 	        {return Struct(i)->mult;}
  Double_t 	L0(Int_t i = 0) 	        {return Struct(i)->L0;}
  Double_t 	bbcX(Int_t i = 0) 	        {return Struct(i)->bbcX;}
  Double_t 	bbcXctbTOFp(Int_t i = 0) 	{return Struct(i)->bbcXctbTOFp;}
  Double_t 	bbcWest(Int_t i = 0) 	        {return Struct(i)->bbcWest;}
  Double_t 	bbcEast(Int_t i = 0) 	        {return Struct(i)->bbcEast;}
  Double_t 	bbcYellowBkg(Int_t i = 0) 	{return Struct(i)->bbcYellowBkg;}
  Double_t 	bbcBlueBkg(Int_t i = 0) 	{return Struct(i)->bbcBlueBkg;}
  Double_t 	pvpdWest(Int_t i = 0) 	        {return Struct(i)->pvpdWest;}
  Double_t 	pvpdEast(Int_t i = 0) 	        {return Struct(i)->pvpdEast;}
  Double_t 	zdcCoin(Int_t i = 0)            {return Nc(zdcX(i), zdcEast(i), zdcWest(i));}
  Double_t 	bbcCoin(Int_t i = 0)            {return Nc(bbcX(i), bbcEast(i), bbcWest(i));}
  void		validityMargin(Double_t margin = 0) {fMargin = margin;}
  Double_t getCTBWest() {return ctbWest();}
  Double_t getCTBEast() {return ctbEast();}
  Double_t getCTBOrTOFp() {return ctbTOFp();}
  Double_t getTOFp() {return tofp();}
  Double_t getZDCWest() {return zdcWest();}
  Double_t getZDCEast() {return zdcEast();}
  Double_t getZDCX() {return zdcX();}
  Double_t getZDCCoin() {return zdcCoin();}
  Double_t getMult() {return mult();}
  Double_t getL0() {return L0();}
  Double_t getBBCX() {return bbcX();}
  Double_t getBBCCoin() {return bbcCoin();}
  Double_t getBBCXCTB() {return bbcXctbTOFp();}
  Double_t getBBCWest() {return bbcWest();}
  Double_t getBBCEast() {return bbcEast();}
  Double_t getBBCYellowBkg() {return bbcYellowBkg();}
  Double_t getBBCBlueBkg() {return bbcBlueBkg();}
  Double_t getPVPDWest() {return pvpdWest();}
  Double_t getPVPDEast() {return pvpdEast();}
  UInt_t   getRichHVStatus() {return St_richvoltagesC::instance()->status();}
  void     setValidityMargin(Double_t margin = 0) {validityMargin(margin);}

  // The following code attempts to correct coincidence rates for accidentals and multiples
  // See STAR Note 528
  static Double_t Nc(Double_t New, Double_t Ne, Double_t Nw, Int_t n_bunches = 111)
  {
    // 111 is a guess using the maximum seen filled bunches in RHIC so far
    // (not always the case, but we don't have access to this number)
    Double_t Nbc = StDetectorDbClock::instance()->CurrentFrequency() * ((Double_t) n_bunches) / 120.;
    return -Nbc * TMath::Log(1. - ((New - (Ne * Nw / Nbc)) / (Nbc + New - Ne - Nw)));
  }
 private:
  static St_trigDetSumsC* fgInstance;
  Double_t	fMargin;

  ClassDefChair(St_trigDetSums, trigDetSums_st )
};
#endif
