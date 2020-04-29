/***********************************************************************
 * Author: brian Feb 6, 1998
 * Description:  Raw data information along with access functions
 ***********************************************************************/
#ifndef ST_TPC_PAD_COORDINATE_HH
#define ST_TPC_PAD_COORDINATE_HH

#include <ostream>

#include "Rtypes.h"
class StTpcPadCoordinate
{
 public:
  StTpcPadCoordinate(const Int_t sector = 0, const Int_t row = 0, const Float_t pad = 0, const Float_t tb = 0) : mSector(sector), mRow(row), mPad(pad), mTimeBucket(tb) {/**/}
  ~StTpcPadCoordinate() {/**/}
  Int_t operator==(const StTpcPadCoordinate &p) const {return (p.mSector == mSector && p.mRow == mRow && p.mPad == mPad && p.mTimeBucket == mTimeBucket);}
  Int_t operator!=(const StTpcPadCoordinate &p) const {return !(*this == p);};
  // access functions
  Int_t sector()           const {return mSector;}
  Int_t row()              const {return mRow;}
  Float_t pad()            const {return mPad;}
  Float_t timeBucket()     const {return mTimeBucket;}
  Int_t sector()                 {return mSector;}
  Int_t row()          		 {return mRow;}
  Float_t pad()          	 {return mPad;}
  Float_t timeBucket()    	 {return mTimeBucket;}

  void setSector(Int_t s)        {mSector = s;}
  void setRow(Int_t r)           {mRow = r;}
  void setPad(Float_t p)           {mPad = p;}
  void setTimeBucket(Float_t t)    {mTimeBucket = t;}

 protected:
  Int_t mSector;
  Int_t mRow;
  Float_t mPad;
  Float_t mTimeBucket;

};
// Non-member
std::ostream &operator<<(std::ostream &, const StTpcPadCoordinate &);

#endif
