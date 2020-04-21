/***************************************************************************
 * Author:  David Hardtke
 *
 * Description: This is the interface between the database and the offline
 *              TPC software.  This classes takes care of the annoying
 *              calls to the root infrastucture, packages and manipulates
 *              the data, and returns the data to the user via simple
 *              interface classes.
 **************************************************************************/
#ifndef ClassStTpcDb
#define ClassStTpcDb

#include "TGeoMatrix.h"

#include "struct_containers.h"

class StTpcDb
{
 public:
  static StTpcDb& instance() {
    static StTpcDb instance;
    return instance;
  }
  // Glob     = Global coordinate
  // Tpc      = Tpc    -"-                survey
  // Half     = Tpc Half west / east -"-  survey
  // SupS     = super sector misalignment(?)
  // SubS[io] = SubSector[io] misalignment
  // SecL     = sector -"- coordinate (y_p, x_p, DriftDistance - z_p);
  // Pad      = Pad -"- (x_p,y_p,z_p) (Sector12 coordinate system)
  // Tpc => Global is mTpc2GlobMatrix
  // Pad => SecL   is internal Flip matrix
  enum ETpcSectorRotationType {kUndefSector     = -2,
                               kFlip            = -1, // Flip * Subs[io] => SupS
                               kSupS2Tpc        = 0, // SupS => Tpc
                               kSupS2Glob       = 1, // SupS => Tpc => Glob;
                               kSubSInner2SupS  = 2, // Subs[io] => SupS
                               kSubSOuter2SupS  = 3, // -"-
                               kSubSInner2Tpc   = 4, // (Subs[io] => SupS) => Tpc
                               kSubSOuter2Tpc   = 5, // -"-
                               kSubSInner2Glob  = 6, // (Subs[io] => SupS => Tpc) => Glob
                               kSubSOuter2Glob  = 7, // -"-
                               kPadInner2SupS   = 8, // (Pad => SecL) => (SubS[io] => SupS)
                               kPadOuter2SupS   = 9, // -"-
                               kPadInner2Tpc    = 10, // (Pad => SecL) => (SubS[io] => SupS => Tpc)
                               kPadOuter2Tpc    = 11, // -"-
                               kPadInner2Glob   = 12, // (Pad => SecL) => (SubS[io] => SupS => Tpc => Glob)
                               kPadOuter2Glob   = 13, // -"-
                               kTotalTpcSectorRotaions = 14
                              };
 private:
  TGeoTranslation*      mSwap[2];       //!
  TGeoHMatrix*          mFlip;          //!
  TGeoHMatrix*          mTpc2GlobMatrix;//!
  TGeoHMatrix*          mHalf[2];       //!
  TGeoHMatrix*          mTpcSectorRotations[24][kTotalTpcSectorRotaions]; //!
  float               mDriftVel[2];   //!
  double              mzGG;           //! Gating Grid z
  static bool         mOldScheme;     //! switch between Old and New alignment scheme
 private:
  StTpcDb();
 public:
  ~StTpcDb();
  St_tpcWirePlanesC*     WirePlaneGeometry() {return St_tpcWirePlanesC::instance();}
  St_tpcDimensionsC*     Dimensions() {return St_tpcDimensionsC::instance();}
  St_tpcElectronicsC*    Electronics() {return St_tpcElectronicsC::instance();}
  St_tpcGlobalPositionC* GlobalPosition() {return St_tpcGlobalPositionC::instance();}
  St_tpcFieldCageC*      FieldCage() {return St_tpcFieldCageC::instance();}
  St_tpcPadResponseC*    PadResponse() {return St_tpcPadResponseC::instance();}
  float                triggerTimeOffset()     {return St_trgTimeOffsetC::instance()->triggerTimeOffset();}
  static bool          IsOldScheme()    {return mOldScheme;}
  void    SetDriftVelocity();
  float DriftVelocity(int sector = 24);
  void SetTpcRotations();
  void SetTpcRotationMatrix(TGeoHMatrix* m, int sector = 0, int k = kSupS2Tpc)
  {
    if (sector == 0)  {if (m) *mTpc2GlobMatrix = *m;}
    else              {if (m) *mTpcSectorRotations[sector - 1][k] = *m;}
  }
  const TGeoHMatrix &Flip()                           const {return *mFlip;}
  const TGeoHMatrix &Tpc2GlobalMatrix()               const {return *mTpc2GlobMatrix;}
  const TGeoHMatrix &TpcRot(int sector, int k)    const {return *mTpcSectorRotations[sector - 1][k];}
  const TGeoHMatrix &SupS2Tpc(int sector = 1)       const {return TpcRot(sector, kSupS2Tpc);}
  const TGeoHMatrix &SupS2Glob(int sector = 1)      const {return TpcRot(sector, kSupS2Glob);}
  const TGeoHMatrix &SubSInner2SupS(int sector = 1) const {return TpcRot(sector, kSubSInner2SupS);}
  const TGeoHMatrix &SubSOuter2SupS(int sector = 1) const {return TpcRot(sector, kSubSOuter2SupS);}
  const TGeoHMatrix &SubSInner2Tpc(int sector = 1)  const {return TpcRot(sector, kSubSInner2Tpc);}
  const TGeoHMatrix &SubSOuter2Tpc(int sector = 1)  const {return TpcRot(sector, kSubSOuter2Tpc);}
  const TGeoHMatrix &SubSInner2Glob(int sector = 1) const {return TpcRot(sector, kSubSInner2Glob);}
  const TGeoHMatrix &SubSOuter2Glob(int sector = 1) const {return TpcRot(sector, kSubSOuter2Glob);}

  const TGeoHMatrix &PadInner2SupS(int sector = 1)  const {return TpcRot(sector, kPadInner2SupS);}
  const TGeoHMatrix &PadOuter2SupS(int sector = 1)  const {return TpcRot(sector, kPadOuter2SupS);}
  const TGeoHMatrix &PadInner2Tpc(int sector = 1)   const {return TpcRot(sector, kPadInner2Tpc);}
  const TGeoHMatrix &PadOuter2Tpc(int sector = 1)   const {return TpcRot(sector, kPadOuter2Tpc);}
  const TGeoHMatrix &PadInner2Glob(int sector = 1)  const {return TpcRot(sector, kPadInner2Glob);}
  const TGeoHMatrix &PadOuter2Glob(int sector = 1)  const {return TpcRot(sector, kPadOuter2Glob);}

  const TGeoHMatrix &SubS2SupS(int sector = 1, int row = 1) const {int k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kSubSInner2SupS : kSubSOuter2SupS; return TpcRot(sector, k);}
  const TGeoHMatrix &SubS2Tpc(int sector = 1, int row = 1)  const {int k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kSubSInner2Tpc : kSubSOuter2Tpc; return TpcRot(sector, k);}
  const TGeoHMatrix &SubS2Glob(int sector = 1, int row = 1) const {int k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kSubSInner2Glob : kSubSOuter2Glob; return TpcRot(sector, k);}

  const TGeoHMatrix &Pad2SupS(int sector = 1, int row = 1)  const {int k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kPadInner2SupS : kPadOuter2SupS; return TpcRot(sector, k);}
  const TGeoHMatrix &Pad2Tpc(int sector = 1, int row = 1)   const {int k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kPadInner2Tpc : kPadOuter2Tpc; return TpcRot(sector, k);}
  const TGeoHMatrix &Pad2Glob(int sector = 1, int row = 1)  const {int k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kPadInner2Glob : kPadOuter2Glob; return TpcRot(sector, k);}
};
#endif
