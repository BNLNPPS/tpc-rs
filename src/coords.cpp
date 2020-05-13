#include <ostream>

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


CoordTransform::CoordTransform()
{
  mTimeBinWidth = 1. / St_tpcElectronicsC::instance()->samplingFrequency();
  mInnerSectorzOffset = St_tpcDimensionsC::instance()->zInnerOffset();
  mOuterSectorzOffset = St_tpcDimensionsC::instance()->zOuterOffset();
}


//      Local Sector Coordnate    <->  Tpc Raw Pad Coordinate
void CoordTransform::operator()(const StTpcLocalSectorCoordinate &a, StTpcPadCoordinate &b, bool useT0, bool useTau)
{
  // useT0 = true for pad and false for cluster, useTau = true for data cluster and  = false for MC
  int sector = a.sector;
  int row    = a.row;

  if (row < 1 || row > St_tpcPadConfigC::instance()->numberOfRows(sector))
    row = rowFromLocalY(a.position.y, a.sector);

  double probablePad = padFromX(a.position.x, a.sector, a.row);
  double                  zoffset = (row > St_tpcPadConfigC::instance()->innerPadRows(sector)) ? mOuterSectorzOffset : mInnerSectorzOffset;
  double t0offset = (useT0 && sector >= 1 && sector <= 24) ? St_tpcPadGainT0BC::instance()->T0(sector, row, tpcrs::irint(probablePad)) : 0;
  t0offset *= mTimeBinWidth;

  if (!useT0 && useTau) // for cluster
    t0offset -= 3.0 * St_tss_tssparC::instance()->tau();   // correct for convolution lagtime

  double t0zoffset = t0offset * StTpcDb::instance().DriftVelocity(sector) * 1e-6;
  double tb = tBFromZ(a.position.z + zoffset - t0zoffset, sector, row, probablePad);
  b = StTpcPadCoordinate{sector, row, probablePad, tb};
}


void CoordTransform::operator()(const StTpcPadCoordinate &a, StTpcLocalSectorCoordinate &b, bool useT0, bool useTau)
{
  // useT0 = true for pad and false for cluster, useTau = true for data cluster and = false for MC
  Coords  tmp{xFromPad(a.sector, a.row, a.pad), yFromRow(a.sector, a.row), 0};
  double zoffset =  (a.row > St_tpcPadConfigC::instance()->innerPadRows(a.sector)) ? mOuterSectorzOffset : mInnerSectorzOffset;
  double t0offset = useT0 ? St_tpcPadGainT0BC::instance()->T0(a.sector, a.row, tpcrs::irint(a.pad)) : 0;
  t0offset *= mTimeBinWidth;

  if (!useT0 && useTau) // for cluster
    t0offset -= 3.0 * St_tss_tssparC::instance()->tau();   // correct for convolution lagtime

  double t0zoffset = t0offset * StTpcDb::instance().DriftVelocity(a.sector) * 1e-6;
  //t0 offset -- DH  27-Mar-00
  double z = zFromTB(a.timeBucket, a.sector, a.row, a.pad) - zoffset + t0zoffset;
  tmp.z = z;
  b = StTpcLocalSectorCoordinate{{tmp.x, tmp.y, tmp.z}, a.sector, a.row};
}


double CoordTransform::padFromX(double x, int sector, int row) const
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  double pitch = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ?
                         St_tpcPadConfigC::instance()->innerSectorPadPitch(sector) :
                         St_tpcPadConfigC::instance()->outerSectorPadPitch(sector);
  // x coordinate in sector 12
  int npads = St_tpcPadConfigC::instance()->numberOfPadsAtRow(sector, row);
  double xL = x;
  int NiRows = St_tpcPadConfigC::instance()->numberOfInnerRows(sector);

  if (NiRows != 13 && row <= NiRows) {
    // iTPC Survey, see  Jim Thomas comments in CoordTransform::xFromPad
    double yRef = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, NiRows) + 0.565;
    double xHit = xL;
    double yHit = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, row) - yRef;
    St_iTPCSurveyC* sur = St_iTPCSurveyC::instance();
    double dx = sur->dx(sector - 1);
    double Xscale = sur->ScaleX(sector - 1);
    double theta  = sur->Angle(sector - 1);
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
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  double pitch = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ?
                   St_tpcPadConfigC::instance()->innerSectorPadPitch(sector) :
                   St_tpcPadConfigC::instance()->outerSectorPadPitch(sector);
  int npads = St_tpcPadConfigC::instance()->numberOfPadsAtRow(sector, row);
  double xPad = -pitch * (pad - (npads + 1.) / 2.);

  int NiRows = St_tpcPadConfigC::instance()->numberOfInnerRows(sector);

  if (NiRows == 13 || row > NiRows) {
    return xPad;
  }

  // iTPC Survey, Jim Thomas correction 08/21/18
  // The change in the yRef comes about because the origin of the coordinate system is 0.565 mm above the center of PadRow 40.
  // The changes for the X coordinates come about because of the reversal of pad counting by the DAQ guys … as you know.
  double yRef = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, NiRows) + 0.565; // Change sign in front of 0.565
  double xL = xPad;  // Eliminate -1 in front of xPad
  double yL = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, row) - yRef;
  St_iTPCSurveyC* sur = St_iTPCSurveyC::instance();
  double dx = sur->dx(sector - 1);
  //double dy = sur->dy(sector-1);
  double Xscale = sur->ScaleX(sector - 1);
  //double Yscale = sur->ScaleY(sector-1);
  double theta  = sur->Angle(sector - 1);
  double xHit = xL * (1. + Xscale) + dx - theta * yL; // Eliminate -1 in front of whole expression and delete ( )
  //double yHit = yL*(1. + Yscale) + dy + theta*xL + yRef;
  return xHit;
}
// Coordinate from Row
//
//Local Transformation...


double CoordTransform::zFromTB(double tb, int sector, int row, int pad) const
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector))
    row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  double trigT0 = StTpcDb::instance().triggerTimeOffset() * 1e6; // units are s
  double elecT0 = St_tpcElectronicsC::instance()->tZero();    // units are us
  double sectT0 = St_tpcPadrowT0C::instance()->T0(sector, row);  // units are us
  double t0 = trigT0 + elecT0 + sectT0;
  int l = sector;

  if ( St_tpcPadConfigC::instance()->IsRowInner(sector, row)) l += 24;

  double tbx = tb + St_tpcSectorT0offsetC::instance()->t0offset(l);

  if (St_tpcRDOT0offsetC::instance()->IsShfited(sector)) {
    tbx += St_tpcRDOT0offsetC::instance()->T0(sector, row, pad);
  }

  double time = t0 + tbx * mTimeBinWidth;
  double z = StTpcDb::instance().DriftVelocity(sector) * 1e-6 * time;
  return z;
}


double CoordTransform::tBFromZ(double z, int sector, int row, int pad) const
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector))
    row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  double trigT0 = StTpcDb::instance().triggerTimeOffset() * 1e6; // units are s
  double elecT0 = St_tpcElectronicsC::instance()->tZero();    // units are us
  double sectT0 = St_tpcPadrowT0C::instance()->T0(sector, row);  // units are us
  double t0 = trigT0 + elecT0 + sectT0;
  double time = z / (StTpcDb::instance().DriftVelocity(sector) * 1e-6);
  int l = sector;

  if ( St_tpcPadConfigC::instance()->IsRowInner(sector, row)) l += 24;

  double tb = (time - t0) / mTimeBinWidth - St_tpcSectorT0offsetC::instance()->t0offset(l);

  if (St_tpcRDOT0offsetC::instance()->IsShfited(sector)) {
    tb -= St_tpcRDOT0offsetC::instance()->T0(sector, row, pad);
  }

  return tb;
}


int CoordTransform::rowFromLocalY(double y, int sector) const
{
  static int Nrows = 0;
  static double* Radii = 0;

  if (Nrows != St_tpcPadConfigC::instance()->padRows(sector)) {
    Nrows = St_tpcPadConfigC::instance()->padRows(sector);

    if (Radii) delete [] Radii;

    Radii = new double[Nrows + 1];

    for (int i = 1; i <= Nrows + 1; i++) {
      if (i == 1) {
        Radii[i - 1] =  (3 * St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, i)
                           - St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, i + 1)) / 2;
      }
      else if (i == Nrows + 1) {
        Radii[i - 1] =  (3 * St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, i - 1)
                           - St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, i - 2)) / 2;
      }
      else {
        Radii[i - 1] = (St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, i - 1) +
                        St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, i)) / 2;
      }
    }
  }

  double* r_ptr = std::lower_bound(Radii, Radii + Nrows + 1, y);
  int row = (r_ptr != Radii + Nrows + 1) && (*r_ptr == y) ? r_ptr - Radii + 1: r_ptr - Radii;

  if (row <= 0) row = 1;
  if (row > Nrows) row = Nrows;

  return row;
}


void  CoordTransform::operator()(const StTpcLocalSectorCoordinate &a, StTpcLocalCoordinate &b)
{
  int row    = a.row;
  int sector = a.sector;

  if (row < 1 || row > St_tpcPadConfigC::instance()->numberOfRows(sector))
    row = rowFromLocalY(a.position.y, a.sector);

  Coords xGG;
  StTpcDb::instance().Pad2Tpc(a.sector, row).LocalToMasterVect(a.position.xyz(), xGG.xyz());
  const double* trans = StTpcDb::instance().Pad2Tpc(sector, row).GetTranslation(); // 4
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  GG2TPC.LocalToMaster(xGG.xyz(), b.position.xyz());
  b.sector = a.sector;
  b.row = row;
}


void  CoordTransform::operator()(const StTpcLocalCoordinate &a, StTpcLocalSectorCoordinate &b)
{
  int row    = a.row;
  int sector = a.sector;

  if (row < 1 || row > St_tpcPadConfigC::instance()->numberOfRows(sector)) {
    Coords xyzS;
    StTpcDb::instance().SupS2Tpc(sector).MasterToLocalVect(a.position.xyz(), xyzS.xyz());
    row = rowFromLocalY(xyzS.x, sector);
  }

  const double* trans = StTpcDb::instance().Pad2Tpc(a.sector, row).GetTranslation(); // 4
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  Coords xGG;
  GG2TPC.MasterToLocal(a.position.xyz(), xGG.xyz());
  StTpcDb::instance().Pad2Tpc(a.sector, row).MasterToLocalVect(xGG.xyz(), b.position.xyz());
  b.sector = a.sector;
  b.row = row;
}
