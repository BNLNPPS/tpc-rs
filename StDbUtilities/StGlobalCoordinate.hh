/***********************************************************************
 * Author:  brian Feb 6, 1998
 *
 * Description:  Raw data information along with access functions
 ***********************************************************************/

#ifndef ST_GLOBAL_COORDINATE_HH
#define ST_GLOBAL_COORDINATE_HH

#include <ostream>

#include "StarClassLibrary/StThreeVector.hh"

class StGlobalCoordinate
{
 public:
  StGlobalCoordinate() {}
  StGlobalCoordinate(const double x, const double y, const double z) : mPosition(x, y, z) { }
  StGlobalCoordinate(const double* x) : mPosition(x) { }
  StGlobalCoordinate(const StThreeVector<double> &x) : mPosition(x) {}
  StGlobalCoordinate(const StThreeVector<float> &x) : mPosition(x.x(), x.y(), x.z()) {}

  ~StGlobalCoordinate() {}
  int operator==(const StGlobalCoordinate &p) const {return p.mPosition == mPosition;}
  int operator!=(const StGlobalCoordinate &p) const {return !(*this == p);}
  // access functions provided by StThreeVector
  const StThreeVector<double> &position() const {return * &mPosition;}
        StThreeVector<double> &position()       {return * &mPosition;}
  void setPosition(const StThreeVector<double> &val) {mPosition = val; }

 protected:
  StThreeVector<double> mPosition;

};
// Non-Member
std::ostream &operator<<(std::ostream &, const StGlobalCoordinate &);
#endif
