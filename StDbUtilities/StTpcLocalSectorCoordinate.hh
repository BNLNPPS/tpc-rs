#ifndef ST_TPC_LOCAL_SECTOR_COORDINATE_HH
#define ST_TPC_LOCAL_SECTOR_COORDINATE_HH

#include <ostream>

#include "StDbUtilities/StTpcCoordinate.h"

class StTpcLocalSectorCoordinate : public StTpcCoordinate
{
 public:
  StTpcLocalSectorCoordinate() : StTpcCoordinate(0, 0, 0, 0, 0) {}
  StTpcLocalSectorCoordinate(double x, double y, double z) :
    StTpcCoordinate(x, y, z, 0, 0) {}
  StTpcLocalSectorCoordinate(const StThreeVector<double> &xyz) :
    StTpcCoordinate(xyz, 0, 0) {}
  StTpcLocalSectorCoordinate(double x, double y, double z, int sector, int row = 0) :
    StTpcCoordinate(x, y, z, sector, row) {}
  StTpcLocalSectorCoordinate(const StThreeVector<double> &xyz, int sector, int row = 0) :
    StTpcCoordinate(xyz, sector, row) {}
  virtual ~StTpcLocalSectorCoordinate() {}
};
std::ostream &operator<<(std::ostream &, const StTpcLocalSectorCoordinate &);
#endif
