/***********************************************************************
 * Author: brian Feb 6, 1998
 *
 * Description:
 *
 * Geometrical transformation Routines for:
 * Raw Pad Coordinate  <-->  Local Coordinate
 *   Local Coordinate  <-->  Global Coordinate
 *
 * These Routines deal positions ONLY!
 ***********************************************************************/
#include "StDbUtilities/StTpcCoordinateTransform.hh"
#include <iostream>
#include "StDetectorDbMaker/St_tpcPadrowT0C.h"
#include "StDetectorDbMaker/St_tpcSectorT0offsetC.h"
#include "StDetectorDbMaker/St_tpcRDOT0offsetC.h"
#include "StDetectorDbMaker/St_tss_tssparC.h"
#include "StDetectorDbMaker/St_tpcPadGainT0BC.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
#include "StDetectorDbMaker/St_tpcPadPlanesC.h"
#include "StDetectorDbMaker/St_iTPCSurveyC.h"
#include "StarClassLibrary/StThreeVector.hh"
#include "tpcrs/logger.h"
#include "math_funcs.h"

static Int_t _debug = 0;

StTpcCoordinateTransform::StTpcCoordinateTransform(StTpcDb* /* globalDbPointer */)
{
  if (St_tpcPadConfigC::instance()
      && StTpcDb::instance()->Electronics()
     ) {
    mTimeBinWidth = 1. / StTpcDb::instance()->Electronics()->samplingFrequency();
    mInnerSectorzOffset = StTpcDb::instance()->Dimensions()->zInnerOffset();
    mOuterSectorzOffset = StTpcDb::instance()->Dimensions()->zOuterOffset();
  }
  else {
    LOG_ERROR << "StTpcDb IS INCOMPLETE! Cannot contstruct Coordinate transformation.\n";
    assert(St_tpcPadConfigC::instance());
    assert(StTpcDb::instance()->Electronics());
  }
}


//      Local Sector Coordnate    <->  Tpc Raw Pad Coordinate
void StTpcCoordinateTransform::operator()(const StTpcLocalSectorCoordinate &a, StTpcPadCoordinate &b, Bool_t useT0, Bool_t useTau)
{
  // useT0 = kTRUE for pad and kFALSE for cluster, useTau = kTRUE for data cluster and  = kFALSE for MC
  Int_t sector = a.fromSector();
  Int_t row    = a.fromRow();

  if (row < 1 || row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row    = rowFromLocal(a);

  Double_t probablePad = padFromLocal(a);
  Double_t                  zoffset = (row > St_tpcPadConfigC::instance()->innerPadRows(sector)) ? mOuterSectorzOffset : mInnerSectorzOffset;
  Double_t t0offset = (useT0 && sector >= 1 && sector <= 24) ? St_tpcPadGainT0BC::instance()->T0(sector, row, tpcrs::irint (probablePad)) : 0;
  t0offset *= mTimeBinWidth;

  if (! useT0 && useTau) // for cluster
    t0offset -= 3.0 * St_tss_tssparC::instance()->tau();   // correct for convolution lagtime

  Double_t t0zoffset = t0offset * StTpcDb::instance()->DriftVelocity(sector) * 1e-6;
  Double_t tb = tBFromZ(a.position().z() + zoffset - t0zoffset, sector, row, probablePad);
  b = StTpcPadCoordinate(sector, row, probablePad, tb);
}


void StTpcCoordinateTransform::operator()(const StTpcPadCoordinate &a,  StTpcLocalSectorCoordinate &b, Bool_t useT0, Bool_t useTau)
{
  // useT0 = kTRUE for pad and kFALSE for cluster, useTau = kTRUE for data cluster and = kFALSE for MC
  StThreeVector<double>  tmp = xyFromRow(a);
  Int_t sector = a.sector();
  Double_t                zoffset =  (a.row() > St_tpcPadConfigC::instance()->innerPadRows(sector)) ? mOuterSectorzOffset : mInnerSectorzOffset;
  Double_t t0offset = useT0 ? St_tpcPadGainT0BC::instance()->T0(a.sector(), a.row(), tpcrs::irint(a.pad())) : 0;
  t0offset *= mTimeBinWidth;

  if (! useT0 && useTau) // for cluster
    t0offset -= 3.0 * St_tss_tssparC::instance()->tau();   // correct for convolution lagtime

  Double_t t0zoffset = t0offset * StTpcDb::instance()->DriftVelocity(a.sector()) * 1e-6;
  //t0 offset -- DH  27-Mar-00
  Double_t z = zFromTB(a.timeBucket(), a.sector(), a.row(), a.pad()) - zoffset + t0zoffset;
  tmp.setZ(z);
  b = StTpcLocalSectorCoordinate(tmp, a.sector(), a.row());
}


Double_t StTpcCoordinateTransform::padFromX(Double_t x, Int_t sector, Int_t row) const
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  Double_t pitch = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ?
                   St_tpcPadConfigC::instance()->innerSectorPadPitch(sector) :
                   St_tpcPadConfigC::instance()->outerSectorPadPitch(sector);
  // x coordinate in sector 12
  Int_t npads = St_tpcPadConfigC::instance()->numberOfPadsAtRow(sector, row);
  Double_t xL = x;
  Int_t NiRows = St_tpcPadConfigC::instance()->numberOfInnerRows(sector);

  if (NiRows != 13 && row <= NiRows) {
    // iTPC Survey, see  Jim Thomas comments in StTpcCoordinateTransform::xFromPad
    Double_t yRef = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, NiRows) + 0.565;
    Double_t xHit = xL;
    Double_t yHit = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, row) - yRef;
    St_iTPCSurveyC* sur = St_iTPCSurveyC::instance();
    Double_t dx = sur->dx(sector - 1);
    //  Double_t dy = sur->dy(sector-1);
    Double_t Xscale = sur->ScaleX(sector - 1);
    //  Double_t Yscale = sur->ScaleY(sector-1);
    Double_t theta  = sur->Angle(sector - 1);
    xL = xHit * (1. - Xscale) - dx + theta * yHit;
    //  Double_t yL = yHit*(1. - Yscale) - dy - theta*xHit + yRef;
  }

  Double_t probablePad = (npads + 1.) / 2. - xL / pitch;

  // CAUTION: pad cannot be <1
  if (probablePad < 0.500001) {
    probablePad = 0.500001;
  }

  if (_debug) {
    LOG_INFO << "StTpcCoordinateTransform::padFromX(" << x << "," << sector << "," << row << "); npads = " << npads << ", pitch = " << pitch
         << "\tprobablePad " << probablePad << '\n';
  }

  return (probablePad);
}


Double_t StTpcCoordinateTransform::xFromPad(Int_t sector, Int_t row, Double_t pad) const      // x coordinate in sector 12
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  Double_t pitch = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ?
                   St_tpcPadConfigC::instance()->innerSectorPadPitch(sector) :
                   St_tpcPadConfigC::instance()->outerSectorPadPitch(sector);
  Int_t npads = St_tpcPadConfigC::instance()->numberOfPadsAtRow(sector, row);
  Double_t xPad = -pitch * (pad - (npads + 1.) / 2.);

  if (_debug) {
    LOG_INFO << "StTpcCoordinateTransform::xFromPad(" << sector << "," << row << "," << pad << "); npads = " << npads << ", pitch = " << pitch
         << "\txPad = " << xPad << '\n';
  }

  Int_t NiRows = St_tpcPadConfigC::instance()->numberOfInnerRows(sector);

  if (NiRows == 13 || row > NiRows) {
    return xPad;
  }

  // iTPC Survey, Jim Thomas correction 08/21/18
  // The change in the yRef comes about because the origin of the coordinate system is 0.565 mm above the center of PadRow 40.
  // The changes for the X coordinates come about because of the reversal of pad counting by the DAQ guys … as you know.
  Double_t yRef = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, NiRows) + 0.565; // Change sign in front of 0.565
  Double_t xL = xPad;  // Eliminate -1 in front of xPad
  Double_t yL = St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, row) - yRef;
  St_iTPCSurveyC* sur = St_iTPCSurveyC::instance();
  Double_t dx = sur->dx(sector - 1);
  //Double_t dy = sur->dy(sector-1);
  Double_t Xscale = sur->ScaleX(sector - 1);
  //Double_t Yscale = sur->ScaleY(sector-1);
  Double_t theta  = sur->Angle(sector - 1);
  Double_t xHit = xL * (1. + Xscale) + dx - theta * yL; // Eliminate -1 in front of whole expression and delete ( )
  //Double_t yHit = yL*(1. + Yscale) + dy + theta*xL + yRef;
  return xHit;
}
// Coordinate from Row
//
//Local Transformation...


Double_t StTpcCoordinateTransform::zFromTB(Double_t tb, Int_t sector, Int_t row, Int_t pad) const
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  Double_t trigT0 = StTpcDb::instance()->triggerTimeOffset() * 1e6;       // units are s
  Double_t elecT0 = StTpcDb::instance()->Electronics()->tZero();          // units are us
  Double_t sectT0 = St_tpcPadrowT0C::instance()->T0(sector, row); // units are us
  Double_t t0 = trigT0 + elecT0 + sectT0;
  Int_t l = sector;

  if ( St_tpcPadConfigC::instance()->IsRowInner(sector, row)) l += 24;

  Double_t tbx = tb + St_tpcSectorT0offsetC::instance()->t0offset(l);

  if (St_tpcRDOT0offsetC::instance()->IsShfited(sector)) {
    tbx += St_tpcRDOT0offsetC::instance()->T0(sector, row, pad);
  }

  Double_t time = t0 + tbx * mTimeBinWidth;
  Double_t z = StTpcDb::instance()->DriftVelocity(sector) * 1e-6 * time;
  return z;
}


Double_t StTpcCoordinateTransform::tBFromZ(Double_t z, Int_t sector, Int_t row, Int_t pad) const
{
  if (row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row = St_tpcPadConfigC::instance()->numberOfRows(sector);

  Double_t trigT0 = StTpcDb::instance()->triggerTimeOffset() * 1e6;       // units are s
  Double_t elecT0 = StTpcDb::instance()->Electronics()->tZero();          // units are us
  Double_t sectT0 = St_tpcPadrowT0C::instance()->T0(sector, row); // units are us
  Double_t t0 = trigT0 + elecT0 + sectT0;
  Double_t time = z / (StTpcDb::instance()->DriftVelocity(sector) * 1e-6);
  Int_t l = sector;

  if ( St_tpcPadConfigC::instance()->IsRowInner(sector, row)) l += 24;

  Double_t tb = (time - t0) / mTimeBinWidth - St_tpcSectorT0offsetC::instance()->t0offset(l);

  if (St_tpcRDOT0offsetC::instance()->IsShfited(sector)) tb -= St_tpcRDOT0offsetC::instance()->T0(sector, row, pad);

  return tb;
}


// FOR SECTOR 12 ONLY!!!! (Local coordinate);
Int_t StTpcCoordinateTransform::rowFromLocalY(Double_t y, Int_t sector)
{
  static Int_t Nrows = 0;
  static Double_t* Radii = 0;
#ifndef __OLD__

  if (Nrows != St_tpcPadConfigC::instance()->padRows(sector)) {
    Nrows = St_tpcPadConfigC::instance()->padRows(sector);

    if (Radii) delete [] Radii;

    Radii = new Double_t[Nrows + 1];

    for (Int_t i = 1; i <= Nrows + 1; i++) {
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
#else

  if (! Nrows) {
    Nrows = St_tpcPadPlanesC::instance()->padRows();
    Radii = new Double_t[Nrows];

    for (Int_t i = 1; i <= Nrows; i++) {
      Radii[i - 1] = St_tpcPadPlanesC::instance()->radialDistanceAtRow(i);
    }
  }

  if (y < Radii[0]) return 1;

  if (y > Radii[Nrows - 1]) return Nrows;

  double* r_ptr = std::lower_bound(Radii, Radii + Nrows, y);
  int row = (r_ptr != Radii + Nrows) && (*r_ptr == y) ? r_ptr - Radii : r_ptr - Radii - 1;

  if (row < Nrows - 1) {
    if (std::abs(Radii[row] - y) > std::abs(Radii[row + 1] - y)) row++;
  }

  row++;
  return row;
#endif
}


void  StTpcCoordinateTransform::operator()(const        StTpcLocalSectorCoordinate &a, StTpcLocalCoordinate &b           )
{
  StThreeVector<double> xGG;
  Int_t row    = a.fromRow();
  Int_t sector = a.fromSector();

  if (row < 1 || row > St_tpcPadConfigC::instance()->numberOfRows(sector)) row    = rowFromLocal(a);

  StTpcDb::instance()->Pad2Tpc(a.sector(), row).LocalToMasterVect(a.position().xyz(), xGG.xyz());
  const Double_t* trans = StTpcDb::instance()->Pad2Tpc(sector, row).GetTranslation(); // 4
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  GG2TPC.LocalToMaster(xGG.xyz(), b.position().xyz());
  b.setSector(a.sector()); b.setRow(row);
}


void  StTpcCoordinateTransform::operator()(const              StTpcLocalCoordinate &a, StTpcLocalSectorCoordinate &b     )
{
  Int_t row    = a.fromRow();
  Int_t sector = a.fromSector();

  if ( ! (row >= 1 && row <= St_tpcPadConfigC::instance()->numberOfRows(sector))) {
    StThreeVector<double> xyzS;
    StTpcDb::instance()->SupS2Tpc(sector).MasterToLocalVect(a.position().xyz(), xyzS.xyz());
    row = rowFromLocalY(xyzS[0], sector);
  }

  const Double_t* trans = StTpcDb::instance()->Pad2Tpc(a.sector(), row).GetTranslation(); // 4
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  StThreeVector<double> xGG;
  GG2TPC.MasterToLocal(a.position().xyz(), xGG.xyz());
  StTpcDb::instance()->Pad2Tpc(a.sector(), row).MasterToLocalVect(xGG.xyz(), b.position().xyz()); b.setSector(a.sector()); b.setRow(row);
}
