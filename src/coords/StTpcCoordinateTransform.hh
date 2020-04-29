/***********************************************************************
 * Author: brian made this on  Feb 6, 1998
 *
 * Description:
 *
 * Geometrical transformation Routines for:
 * Raw Pad Coordinate  <-->  Local Coordinate
 *   Local Coordinate  <-->  Global Coordinate
 *
 * These Routines deal positions ONLY!
 ***********************************************************************/
#ifndef ST_COORDINATE_TRANSFORM_HH
#define ST_COORDINATE_TRANSFORM_HH

#include <stdlib.h>
// SCL
#include "particles/SystemOfUnits.h"
#include "particles/StThreeVector.hh"

//#define DEBUG_TPC 0
//#define idb if(DEBUG_TPC) LOG_INFO
#include "coords/StTpcPadCoordinate.hh"
#include "coords/StTpcLocalCoordinate.hh"
#include "coords/StTpcLocalSectorCoordinate.hh"
#include "coords/StGlobalCoordinate.hh"

#include "coords/StTpcLocalDirection.hh"
#include "coords/StTpcLocalSectorDirection.hh"
#include "coords/StGlobalDirection.hh"

#include "tpc_db.h"

// pad              => sector12       =>   subsector => sector => tpc      => global
// TpcPadCoordinate => TpcSectL => TpcSectLAligned => TpcLocal => Global
class StTpcCoordinateTransform
{
 public:
  StTpcCoordinateTransform(StTpcDb* globalDbPointer = 0);
  ~StTpcCoordinateTransform() {}
  //      Raw Data          <--> Tpc Local Sector Coordinates
  void  operator()(const StTpcLocalSectorCoordinate &a, StTpcPadCoordinate &b, Bool_t useT0 = kFALSE, Bool_t useTau = kTRUE);
  void  operator()(const StTpcPadCoordinate &a, StTpcLocalSectorCoordinate &b, Bool_t useT0 = kFALSE, Bool_t useTau = kTRUE);
  //      Raw Data          <--> Tpc Local  Coordinates
  void  operator()(const StTpcLocalCoordinate &a, StTpcPadCoordinate &b, Bool_t useT0 = kFALSE, Bool_t useTau = kTRUE)
  {StTpcLocalSectorCoordinate c; this->operator()(a, c); this->operator()(c, b, useT0, useTau);}
  void  operator()(const StTpcPadCoordinate &a, StTpcLocalCoordinate &b, Bool_t useT0 = kFALSE, Bool_t useTau = kTRUE)
  {StTpcLocalSectorCoordinate c; this->operator()(a, c, useT0, useTau); this->operator()(c, b);}
  // Tpc Local Sector <--> TPC Local
  void  operator()(const        StTpcLocalSectorCoordinate &a, StTpcLocalCoordinate &b           );
  //   { Double_t xyzS[3] = {a.position().x(), a.position().y(), 0.};
  //     StTpcDb::instance()->Pad2Tpc(a.sector(),a.row()).LocalToMaster(xyzS,b.position().xyz());
  //     b.position().setZ(b.position().z() + a.position().z());
  //     b.setSector(a.sector()); b.setRow(a.row());}
  void  operator()(const        StTpcLocalSectorDirection  &a, StTpcLocalDirection               &b)
  {StTpcDb::instance()->Pad2Tpc(a.sector(), a.row()).LocalToMasterVect(a.position().xyz(), b.position().xyz()); b.setSector(a.sector()); b.setRow(a.row());}
  void  operator()(const        StTpcLocalSectorCoordinate &a, StGlobalCoordinate &b)
  { StTpcLocalCoordinate c;    this->operator()(a, c);    this->operator()(c, b);  }
  void  operator()(const  StTpcLocalSectorDirection &a, StGlobalDirection &b)
  {StTpcDb::instance()->Pad2Glob(a.sector(), a.row()).LocalToMasterVect(a.position().xyz(), b.position().xyz());}
  void  operator()(const              StTpcLocalCoordinate &a, StTpcLocalSectorCoordinate &b     );
  //   { Double_t xyzS[3] = {a.position().x(), a.position().y(), 0.};
  //     StTpcDb::instance()->Pad2Tpc(a.sector(),a.row()).MasterToLocal(xyzS,b.position().xyz()); b.setSector(a.sector()); b.setRow(a.row());
  //     b.position().setZ(b.position().z() + a.position().z());
  //   }
  void  operator()(const               StTpcLocalDirection &a, StTpcLocalSectorDirection &b      )
  {StTpcDb::instance()->Pad2Tpc(a.sector(), a.row()).MasterToLocalVect(a.position().xyz(), b.position().xyz()); b.setSector(a.sector()); b.setRow(a.row());}
  void  operator()(const               StGlobalCoordinate &a, StTpcLocalSectorCoordinate &b, Int_t sector, Int_t row)
  {
    StTpcLocalCoordinate c;
    this->operator()(a, c, sector, row);
    this->operator()(c, b);
  }
  void  operator()(const  StGlobalDirection &a, StTpcLocalSectorDirection &b, Int_t sector, Int_t row)
  {StTpcDb::instance()->Pad2Glob(sector, row).MasterToLocalVect(a.position().xyz(), b.position().xyz()); b.setSector(sector); b.setRow(row);}
  // Internal TpcCoordinate <-->  Global Coordinate
  void  operator()(const StTpcLocalCoordinate &a, StGlobalCoordinate &b)
  {StTpcDb::instance()->Tpc2GlobalMatrix().LocalToMaster(a.position().xyz(), b.position().xyz());}
  void  operator()(const StGlobalCoordinate &a, StTpcLocalCoordinate &b, Int_t sector, Int_t row)
  {StTpcDb::instance()->Tpc2GlobalMatrix().MasterToLocal(a.position().xyz(), b.position().xyz()); b.setSector(sector); b.setRow(row);}
  void  operator()(const StTpcLocalDirection &a, StGlobalDirection &b)
  {StTpcDb::instance()->Tpc2GlobalMatrix().LocalToMasterVect(a.position().xyz(), b.position().xyz());}
  void  operator()(const StGlobalDirection &a, StTpcLocalDirection &b, Int_t sector, Int_t row)
  {StTpcDb::instance()->Tpc2GlobalMatrix().MasterToLocalVect(a.position().xyz(), b.position().xyz()); b.setSector(sector); b.setRow(row);}
  //      Raw Data          <-->  Global Coordinate
  void  operator()(const StTpcPadCoordinate &a, StGlobalCoordinate &b, Bool_t useT0 = kFALSE, Bool_t useTau = kTRUE)
  {StTpcLocalCoordinate c; this->operator()(a, c, useT0, useTau); this->operator()(c, b);}
  void  operator()(const StGlobalCoordinate &a, StTpcPadCoordinate &b, Int_t sector, Int_t row, Bool_t useT0 = kFALSE, Bool_t useTau = kTRUE)
  {StTpcLocalCoordinate c; this->operator()(a, c, sector, row); this->operator()(c, b, useT0, useTau);}
  Double_t   tBFromZ(Double_t z, Int_t sector, Int_t row, Int_t pad = 0) const;
  Double_t  zFromTB(Double_t tb, Int_t sector, Int_t row, Int_t pad = 0) const;
  // Transformation Routines!!
  // Raw Data (pad row timebin or drift L From tpc local sector Coordinates
  static Int_t       rowFromLocalY(Double_t y, Int_t sector);
  static Int_t       rowFromLocal(const StThreeVector<Double_t> &a, Int_t sector)            {return rowFromLocalY(a.y(), sector);}
  Double_t    padFromLocal(const StThreeVector<Double_t> &a, Int_t sector, Int_t row)  const {return padFromX(a.x(), sector, row);}
  Double_t    padFromX(Double_t x, Int_t sector, Int_t row)                        const;
  Int_t       rowFromLocal(const StTpcLocalSectorCoordinate &a)      const {return rowFromLocal(a.position(), a.sector());}
  Double_t    padFromLocal(const StTpcLocalSectorCoordinate &a)      const {return padFromLocal(a.position(), a.sector(), a.row());}
  // tpc local sector Coordinates from Raw Data
  StThreeVector<Double_t> xyFromRow(const StTpcPadCoordinate &a) {return StThreeVector<Double_t> (xFromPad(a.sector(), a.row(), a.pad()), yFromRow(a.sector(), a.row()), 0);}
  Double_t                yFromRow(Int_t sector, Int_t row)                        const {return (St_tpcPadConfigC::instance()->radialDistanceAtRow(sector, row));}
  Double_t                xFromPad(Int_t sector, Int_t row, Double_t pad)          const;
  // sector from Tpc local coordinates
  Int_t sectorFromCoordinate(const StThreeVector<double> &a) const
  {
    Double_t angle = 180 / M_PI * std::atan2(a.y(), a.x());

    if (angle < 0) angle += 360;

    Int_t sectorNumber = (int)( (angle + 15) / 30);

    if (a.z() > 0) {sectorNumber = 15 - sectorNumber; if (sectorNumber > 12)sectorNumber -= 12;}
    else       {sectorNumber += 9;              if (sectorNumber <= 12)sectorNumber += 12;}

    return sectorNumber;
  }
 private:
  Double_t    mTimeBinWidth;
  Double_t    mInnerSectorzOffset;
  Double_t    mOuterSectorzOffset;
  Int_t       mNoOfInnerRows;
  Int_t       mNoOfRows;

};

#endif



