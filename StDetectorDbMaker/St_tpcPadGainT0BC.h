#ifndef St_tpcPadGainT0BC_h
#define St_tpcPadGainT0BC_h
#include "TObject.h"
#include "StDetectorDbMaker/St_tpcPadGainT0C.h"
#include "StDetectorDbMaker/St_itpcPadGainT0C.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
class St_tpcPadGainT0BC : public TObject
{
 public:
  static St_tpcPadGainT0BC* 	instance();
  float 	Gain(int sector, int row, int pad) const;
  float 	  T0(int sector, int row, int pad) const;
  bool    livePadrow(int sector, int row) const;
 protected:
  St_tpcPadGainT0BC() {}
  virtual ~St_tpcPadGainT0BC() {fgInstance = 0;}
 private:
  static St_tpcPadGainT0BC* fgInstance;
};
#endif
