#ifndef ST_GLOBAL_DIRECTION_HH
#define ST_GLOBAL_DIRECTION_HH

#include <ostream>

#include "StDbUtilities/StGlobalCoordinate.hh"
class StGlobalDirection : public StGlobalCoordinate
{
 public:
  StGlobalDirection() : StGlobalCoordinate() {}
  StGlobalDirection(const double x, const double y, const double z) :
    StGlobalCoordinate(x, y, z) {}
  StGlobalDirection(const StThreeVector<double> &xyz) : StGlobalCoordinate(xyz) {}
  StGlobalDirection(const StThreeVector<float> &xyz) :  StGlobalCoordinate(xyz) {}
  ~StGlobalDirection() {};
};
// Non-Member
std::ostream &operator<<(std::ostream &, const StGlobalDirection &);
#endif
