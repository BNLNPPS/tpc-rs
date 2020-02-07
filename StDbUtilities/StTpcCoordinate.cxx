/***********************************************************************
 * $Id: StTpcCoordinate.cxx,v 1.2 2011/01/18 14:34:28 fisyak Exp $
 ***********************************************************************/
#include "StDbUtilities/StTpcLocalDirection.hh"
#include "StDbUtilities/StTpcLocalCoordinate.hh"
#include "StDbUtilities/StTpcLocalSectorAlignedDirection.hh"
#include "StDbUtilities/StTpcLocalSectorAlignedCoordinate.hh"
#include "StDbUtilities/StTpcLocalSectorDirection.hh"
#include "StDbUtilities/StTpcLocalSectorCoordinate.hh"


StTpcLocalCoordinate::StTpcLocalCoordinate(double x, double y, double z) : StTpcCoordinate(x, y, z, 0, 0) {}


StTpcLocalCoordinate::StTpcLocalCoordinate(const StThreeVector<double> &xyz) : StTpcCoordinate(xyz, 0, 0) {}


#define OS "( (" <<  a.position().x() << ", " \
    << a.position().y() << ", " \
    << a.position().z() << ") " \
    << ", " << a.fromSector() << "," << a.fromRow() << " )"
// Non-member Functions
ostream &operator<<(ostream &os, const StTpcCoordinate &a)
{
  return os << OS;
}


ostream &operator<<(ostream &os, const StTpcLocalDirection &a)
{
  return os << "TPC_Local Direction( (" << OS;
}


ostream &operator<<(ostream &os,
                    const StTpcLocalCoordinate &a)
{
  return os << "TPC_Local( (" << OS;
}


ostream &operator<<(ostream &os, const StTpcLocalSectorCoordinate &a)
{
  return os << "TPC_Local_Sector( (" << OS;
}
ostream &operator<<(ostream &os, const StTpcLocalSectorDirection &a)
{
  return os << "TPC_Local_Sector Direction( (" << OS;
}
#undef OS
