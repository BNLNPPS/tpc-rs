#ifndef ST_TPC_LOCAL_COORDINATE_HH
#define ST_TPC_LOCAL_COORDINATE_HH

#include <ostream>

#include "coords/StTpcCoordinate.h"

class StTpcLocalCoordinate : public StTpcCoordinate
{
 public:
  StTpcLocalCoordinate() :  StTpcCoordinate(0, 0, 0, 0, 0) {}
  StTpcLocalCoordinate(double x, double y, double z);
  StTpcLocalCoordinate(double x, double y, double z, int sector, int row) :
    StTpcCoordinate(x, y, z, sector, row) {}
  StTpcLocalCoordinate(const StThreeVector<double> &xyz);
  StTpcLocalCoordinate(const StThreeVector<double> &xyz, int sector, int row) :
    StTpcCoordinate(xyz, sector, row) {}
  ~StTpcLocalCoordinate() {}
};
std::ostream &operator<<(std::ostream &, const StTpcLocalCoordinate &);
#endif
