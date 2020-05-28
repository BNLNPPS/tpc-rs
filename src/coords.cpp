#include <ostream>

#include "tpcrs/configurator.h"
#include "coords.h"
#include "math_funcs.h"

double mag(const Coords& c)
{
  return std::sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
}

double perp(const Coords& c)
{
  return std::sqrt(c.x*c.x + c.y*c.y);
}

Coords unit(const Coords& c)
{
  double m = mag(c);
  return m ? Coords{c.x/m, c.y/m, c.z/m} : Coords{};
}

Coords operator-(const Coords &c1, const Coords &c2)
{
  return Coords{c1.x - c2.x, c1.y - c2.y, c1.z - c2.z};
}

double operator*(const Coords &c1, const Coords &c2)
{
  return c1.x*c2.x + c1.y*c2.y + c1.z*c2.z;
}

Coords operator/(const Coords &c, double v)
{
  return Coords{c.x/v, c.y/v, c.z/v};
}

std::ostream &operator<<(std::ostream &os, const StTpcPadCoordinate &a)
{
  return os << "(sector= " << a.sector
         << ", row= "    << a.row
         << ", pad= "    << a.pad
         << ", tbuck= "  << a.timeBucket << ")";
}

std::ostream &operator<<(std::ostream &os, const StGlobalCoordinate &a)
{
  return os << "GC ( "
         << a.position.x << ", "
         << a.position.y << ", "
         << a.position.z << ")";
}

std::ostream &operator<<(std::ostream &os, const StGlobalDirection &a)
{
  return os << "GD ( "
         << a.position.x << ", "
         << a.position.y << ", "
         << a.position.z << ")";
}

#define OS "( (" <<  a.position.x << ", " \
    << a.position.y << ", " \
    << a.position.z << ") " \
    << ", " << a.sector << "," << a.row << " )"

std::ostream &operator<<(std::ostream &os, const StTpcCoordinate &a)
{
  return os << OS;
}

std::ostream &operator<<(std::ostream &os, const StTpcLocalCoordinate &a)
{
  return os << "TPC_Local( (" << OS;
}

std::ostream &operator<<(std::ostream &os, const StTpcLocalDirection &a)
{
  return os << "TPC_Local Direction( (" << OS;
}

std::ostream &operator<<(std::ostream &os, const StTpcLocalSectorCoordinate &a)
{
  return os << "TPC_Local_Sector( (" << OS;
}

std::ostream &operator<<(std::ostream &os, const StTpcLocalSectorDirection &a)
{
  return os << "TPC_Local_Sector Direction( (" << OS;
}

using tpcrs::Cfg;

CoordTransform::CoordTransform(const tpcrs::Configurator& cfg) :
  cfg_(cfg)
{
  mTimeBinWidth = 1e6 / cfg_.S<starClockOnl>().frequency;
  mInnerSectorzOffset = cfg_.S<tpcEffectiveGeom>().z_inner_offset;
  mOuterSectorzOffset = cfg_.S<tpcEffectiveGeom>().z_outer_offset;
}


// Local Sector Coordnate <-> Tpc Raw Pad Coordinate
void CoordTransform::local_sector_to_hardware(const StTpcLocalSectorCoordinate &a, StTpcPadCoordinate &b, bool useT0, bool useTau) const
{
  // useT0 = true for pad and false for cluster, useTau = true for data cluster and  = false for MC
  int row = a.row;

  if (row < 1 || row > cfg_.C<St_tpcPadConfigC>().numberOfRows(a.sector))
    row = rowFromLocalY(a.position.y, a.sector);

  double probablePad = padFromX(a.position.x, a.sector, a.row);
  double zoffset = (row > cfg_.C<St_tpcPadConfigC>().innerPadRows(a.sector)) ? mOuterSectorzOffset : mInnerSectorzOffset;
  double t0offset = (useT0 && a.sector >= 1 && a.sector <= 24) ? cfg_.C<St_tpcPadGainT0BC>().T0(a.sector, row, tpcrs::irint(probablePad)) : 0;
  t0offset *= mTimeBinWidth;

  if (!useT0 && useTau) // for cluster
    t0offset -= 3.0 * cfg_.S<tss_tsspar>().tau;   // correct for convolution lagtime

  double t0zoffset = t0offset * tpcrs::DriftVelocity(a.sector, cfg_) * 1e-6;
  double tb = tBFromZ(a.position.z + zoffset - t0zoffset, a.sector, row, probablePad);
  b = StTpcPadCoordinate{a.sector, row, probablePad, tb};
}


void CoordTransform::hardware_to_local_sector(const StTpcPadCoordinate &a, StTpcLocalSectorCoordinate &b, bool useT0, bool useTau) const
{
  // useT0 = true for pad and false for cluster, useTau = true for data cluster and = false for MC
  Coords  tmp{xFromPad(a.sector, a.row, a.pad), yFromRow(a.sector, a.row), 0};
  double zoffset =  (a.row > cfg_.C<St_tpcPadConfigC>().innerPadRows(a.sector)) ? mOuterSectorzOffset : mInnerSectorzOffset;
  double t0offset = useT0 ? cfg_.C<St_tpcPadGainT0BC>().T0(a.sector, a.row, tpcrs::irint(a.pad)) : 0;
  t0offset *= mTimeBinWidth;

  if (!useT0 && useTau) // for cluster
    t0offset -= 3.0 * cfg_.S<tss_tsspar>().tau;   // correct for convolution lagtime

  double t0zoffset = t0offset * tpcrs::DriftVelocity(a.sector, cfg_) * 1e-6;
  //t0 offset -- DH  27-Mar-00
  double z = zFromTB(a.timeBucket, a.sector, a.row, a.pad) - zoffset + t0zoffset;
  tmp.z = z;
  b = StTpcLocalSectorCoordinate{{tmp.x, tmp.y, tmp.z}, a.sector, a.row};
}


double CoordTransform::padFromX(double x, int sector, int row) const
{
  if (row > cfg_.C<St_tpcPadConfigC>().numberOfRows(sector)) row = cfg_.C<St_tpcPadConfigC>().numberOfRows(sector);

  double pitch = (row <= cfg_.C<St_tpcPadConfigC>().innerPadRows(sector)) ?
                         cfg_.C<St_tpcPadConfigC>().innerSectorPadPitch(sector) :
                         cfg_.C<St_tpcPadConfigC>().outerSectorPadPitch(sector);
  // x coordinate in sector 12
  int npads = cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(sector, row);
  double xL = x;
  int NiRows = cfg_.C<St_tpcPadConfigC>().numberOfInnerRows(sector);

  if (NiRows != 13 && row <= NiRows) {
    // iTPC Survey, see  Jim Thomas comments in CoordTransform::xFromPad
    double yRef = cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, NiRows) + 0.565;
    double xHit = xL;
    double yHit = cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, row) - yRef;
    const iTPCSurvey& sur = cfg_.S<iTPCSurvey>(sector - 1);
    double dx = sur.dx;
    double Xscale = sur.ScaleX;
    double theta  = sur.Angle;
    xL = xHit * (1. - Xscale) - dx + theta * yHit;
  }

  double probablePad = (npads + 1.) / 2. - xL / pitch;

  // CAUTION: pad cannot be <1
  if (probablePad < 0.500001) {
    probablePad = 0.500001;
  }

  return (probablePad);
}


double CoordTransform::xFromPad(int sector, int row, double pad) const      // x coordinate in sector 12
{
  if (row > cfg_.C<St_tpcPadConfigC>().numberOfRows(sector)) row = cfg_.C<St_tpcPadConfigC>().numberOfRows(sector);

  double pitch = (row <= cfg_.C<St_tpcPadConfigC>().innerPadRows(sector)) ?
                   cfg_.C<St_tpcPadConfigC>().innerSectorPadPitch(sector) :
                   cfg_.C<St_tpcPadConfigC>().outerSectorPadPitch(sector);
  int npads = cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(sector, row);
  double xPad = -pitch * (pad - (npads + 1.) / 2.);

  int NiRows = cfg_.C<St_tpcPadConfigC>().numberOfInnerRows(sector);

  if (NiRows == 13 || row > NiRows) {
    return xPad;
  }

  // iTPC Survey, Jim Thomas correction 08/21/18
  // The change in the yRef comes about because the origin of the coordinate system is 0.565 mm above the center of PadRow 40.
  // The changes for the X coordinates come about because of the reversal of pad counting by the DAQ guys … as you know.
  double yRef = cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, NiRows) + 0.565; // Change sign in front of 0.565
  double xL = xPad;  // Eliminate -1 in front of xPad
  double yL = cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, row) - yRef;
  const iTPCSurvey& sur = cfg_.S<iTPCSurvey>(sector-1);
  double dx = sur.dx;
  double Xscale = sur.ScaleX;
  double theta  = sur.Angle;
  double xHit = xL * (1. + Xscale) + dx - theta * yL; // Eliminate -1 in front of whole expression and delete ( )
  //double yHit = yL*(1. + Yscale) + dy + theta*xL + yRef;
  return xHit;
}
// Coordinate from Row
//
//Local Transformation...


double CoordTransform::zFromTB(double tb, int sector, int row, int pad) const
{
  if (row > cfg_.C<St_tpcPadConfigC>().numberOfRows(sector))
    row = cfg_.C<St_tpcPadConfigC>().numberOfRows(sector);

  double trigT0 = StTpcDb::instance().triggerTimeOffset() * 1e6; // units are s
  double elecT0 = cfg_.S<tpcElectronics>().tZero;    // units are us
  double sectT0 = cfg_.S<tpcPadrowT0>(sector-1).T0[row-1];  // units are us
  double t0 = trigT0 + elecT0 + sectT0;
  int l = sector;

  if ( cfg_.C<St_tpcPadConfigC>().IsRowInner(sector, row)) l += 24;

  double tbx = tb + cfg_.C<St_tpcSectorT0offsetC>().t0offset(l);

  if (cfg_.C<St_tpcRDOT0offsetC>().IsShfited(sector)) {
    tbx += cfg_.C<St_tpcRDOT0offsetC>().T0(sector, row, pad);
  }

  double time = t0 + tbx * mTimeBinWidth;
  double z = tpcrs::DriftVelocity(sector, cfg_) * 1e-6 * time;
  return z;
}


double CoordTransform::tBFromZ(double z, int sector, int row, int pad) const
{
  if (row > cfg_.C<St_tpcPadConfigC>().numberOfRows(sector))
    row = cfg_.C<St_tpcPadConfigC>().numberOfRows(sector);

  double trigT0 = StTpcDb::instance().triggerTimeOffset() * 1e6; // units are s
  double elecT0 = cfg_.S<tpcElectronics>().tZero;    // units are us
  double sectT0 = cfg_.S<tpcPadrowT0>(sector-1).T0[row-1];  // units are us
  double t0 = trigT0 + elecT0 + sectT0;
  double time = z / (tpcrs::DriftVelocity(sector, cfg_) * 1e-6);
  int l = sector;

  if ( cfg_.C<St_tpcPadConfigC>().IsRowInner(sector, row)) l += 24;

  double tb = (time - t0) / mTimeBinWidth - cfg_.C<St_tpcSectorT0offsetC>().t0offset(l);

  if (cfg_.C<St_tpcRDOT0offsetC>().IsShfited(sector)) {
    tb -= cfg_.C<St_tpcRDOT0offsetC>().T0(sector, row, pad);
  }

  return tb;
}


int CoordTransform::rowFromLocalY(double y, int sector) const
{
  static int Nrows = 0;
  static double* Radii = 0;

  if (Nrows != cfg_.C<St_tpcPadConfigC>().padRows(sector)) {
    Nrows = cfg_.C<St_tpcPadConfigC>().padRows(sector);

    if (Radii) delete [] Radii;

    Radii = new double[Nrows + 1];

    for (int i = 1; i <= Nrows + 1; i++) {
      if (i == 1) {
        Radii[i - 1] =  (3 * cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, i)
                           - cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, i + 1)) / 2;
      }
      else if (i == Nrows + 1) {
        Radii[i - 1] =  (3 * cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, i - 1)
                           - cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, i - 2)) / 2;
      }
      else {
        Radii[i - 1] = (cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, i - 1) +
                        cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(sector, i)) / 2;
      }
    }
  }

  // Based on the row position (i.e. "radius") perform a binary search for the
  // row corresponding to the value of y
  double* r_ptr = std::lower_bound(Radii, Radii + Nrows + 1, y);
  int row = (r_ptr != Radii + Nrows + 1) && (*r_ptr == y) ? r_ptr - Radii + 1: r_ptr - Radii;

  if (row <= 0) row = 1;
  if (row > Nrows) row = Nrows;

  return row;
}


void CoordTransform::local_sector_to_local(const StTpcLocalSectorCoordinate &a, StTpcLocalCoordinate &b) const
{
  int row = a.row;

  if (row < 1 || row > cfg_.C<St_tpcPadConfigC>().numberOfRows(a.sector))
    row = rowFromLocalY(a.position.y, a.sector);

  Coords xGG;
  StTpcDb::instance().Pad2Tpc(a.sector, row).LocalToMasterVect(a.position.xyz(), xGG.xyz());
  const double* trans = StTpcDb::instance().Pad2Tpc(a.sector, row).GetTranslation();
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  GG2TPC.LocalToMaster(xGG.xyz(), b.position.xyz());

  b.row = row;
  b.sector = a.sector;
}


void CoordTransform::local_to_local_sector(const StTpcLocalCoordinate &a, StTpcLocalSectorCoordinate &b) const
{
  int row = a.row;

  if (row < 1 || row > cfg_.C<St_tpcPadConfigC>().numberOfRows(a.sector)) {
    Coords xyzS;
    StTpcDb::instance().SupS2Tpc(a.sector).MasterToLocalVect(a.position.xyz(), xyzS.xyz());
    row = rowFromLocalY(xyzS.x, a.sector);
  }

  const double* trans = StTpcDb::instance().Pad2Tpc(a.sector, row).GetTranslation();
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  Coords xGG;
  GG2TPC.MasterToLocal(a.position.xyz(), xGG.xyz());
  StTpcDb::instance().Pad2Tpc(a.sector, row).MasterToLocalVect(xGG.xyz(), b.position.xyz());

  b.row = row;
  b.sector = a.sector;
}
