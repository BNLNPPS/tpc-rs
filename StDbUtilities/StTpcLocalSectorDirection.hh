#ifndef ST_TPC_LOCAL_SECTOR_DIRECTION_HH
#define ST_TPC_LOCAL_SECTOR_DIRECTION_HH

#include <ostream>

#include "StDbUtilities/StTpcCoordinate.h"

class StTpcLocalSectorDirection : public StTpcCoordinate
{
 public:
  StTpcLocalSectorDirection() : StTpcCoordinate(0, 0, 0, 0, 0) {}
  StTpcLocalSectorDirection(double x, double y, double z) :
    StTpcCoordinate(x, y, z, 0, 0) {}
  StTpcLocalSectorDirection(const StThreeVector<double> &xyz) :
    StTpcCoordinate(xyz, 0, 0) {}
  StTpcLocalSectorDirection(double x, double y, double z, int sector, int row = 0) :
    StTpcCoordinate(x, y, z, sector, row) {}
  StTpcLocalSectorDirection(const StThreeVector<double> &xyz, int sector, int row = 0) :
    StTpcCoordinate(xyz, sector, row) {}
  ~StTpcLocalSectorDirection() {}
};
std::ostream &operator<<(std::ostream &, const StTpcLocalSectorDirection &);
#endif
