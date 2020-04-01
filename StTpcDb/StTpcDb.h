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

class StMagUtilities;
#include "St_base/StMessMgr.h"
#include "StEvent/StEnumerations.h"
#include "StDetectorDbMaker/St_tpcPadPlanesC.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
#include "StDetectorDbMaker/St_tpcWirePlanesC.h"
#include "StDetectorDbMaker/St_tpcDimensionsC.h"
#include "StDetectorDbMaker/St_tpcElectronicsC.h"
#include "StDetectorDbMaker/St_tpcSlowControlSimC.h"
#include "StDetectorDbMaker/St_tpcGlobalPositionC.h"
#include "StDetectorDbMaker/St_tpcSectorPositionC.h"
#include "StDetectorDbMaker/St_tpcFieldCageC.h"
#include "StDetectorDbMaker/St_tpcPedestalC.h"
#include "StDetectorDbMaker/St_tpcPadResponseC.h"
#include "StDetectorDbMaker/St_tpcPadGainT0BC.h"
#include "StDetectorDbMaker/St_trgTimeOffsetC.h"
#include "TGeoMatrix.h"
#include "TString.h"
class StTpcDb;
// Global pointers:
R__EXTERN StTpcDb* gStTpcDb;
class StTpcDb
{
 public:
  static StTpcDb* instance() {if (! gStTpcDb) new StTpcDb(); return gStTpcDb;}
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
  StMagUtilities*       mExB;           //!
  TGeoTranslation*      mSwap[2];       //!
  TGeoHMatrix*          mFlip;          //!
  TGeoHMatrix*          mTpc2GlobMatrix;//!
  TGeoHMatrix*          mHalf[2];       //!
  TGeoHMatrix*          mTpcSectorRotations[24][kTotalTpcSectorRotaions]; //!
  Float_t               mDriftVel[2];   //!
  UInt_t                mUc;            //! time for which above mDriftVel have been calculated
  Double_t              mzGG;           //! Gating Grid z
  static Bool_t         mOldScheme;     //! switch between Old and New alignment scheme
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
  Float_t                triggerTimeOffset()     {return St_trgTimeOffsetC::instance()->triggerTimeOffset();}
  static Bool_t          IsOldScheme()    {return mOldScheme;}
  void    SetDriftVelocity();
  Float_t DriftVelocity(Int_t sector = 24);
  void SetTpcRotations();
  void SetTpcRotationMatrix(TGeoHMatrix* m, Int_t sector = 0, Int_t k = kSupS2Tpc)
  {
    if (sector == 0)  {if (m) *mTpc2GlobMatrix = *m;}
    else              {if (m) *mTpcSectorRotations[sector - 1][k] = *m;}
  }
  const TGeoHMatrix &Flip()                           const {return *mFlip;}
  const TGeoHMatrix &Tpc2GlobalMatrix()               const {return *mTpc2GlobMatrix;}
  const TGeoHMatrix &TpcRot(Int_t sector, Int_t k)    const {return *mTpcSectorRotations[sector - 1][k];}
  const TGeoHMatrix &SupS2Tpc(Int_t sector = 1)       const {return TpcRot(sector, kSupS2Tpc);}
  const TGeoHMatrix &SupS2Glob(Int_t sector = 1)      const {return TpcRot(sector, kSupS2Glob);}
  const TGeoHMatrix &SubSInner2SupS(Int_t sector = 1) const {return TpcRot(sector, kSubSInner2SupS);}
  const TGeoHMatrix &SubSOuter2SupS(Int_t sector = 1) const {return TpcRot(sector, kSubSOuter2SupS);}
  const TGeoHMatrix &SubSInner2Tpc(Int_t sector = 1)  const {return TpcRot(sector, kSubSInner2Tpc);}
  const TGeoHMatrix &SubSOuter2Tpc(Int_t sector = 1)  const {return TpcRot(sector, kSubSOuter2Tpc);}
  const TGeoHMatrix &SubSInner2Glob(Int_t sector = 1) const {return TpcRot(sector, kSubSInner2Glob);}
  const TGeoHMatrix &SubSOuter2Glob(Int_t sector = 1) const {return TpcRot(sector, kSubSOuter2Glob);}

  const TGeoHMatrix &PadInner2SupS(Int_t sector = 1)  const {return TpcRot(sector, kPadInner2SupS);}
  const TGeoHMatrix &PadOuter2SupS(Int_t sector = 1)  const {return TpcRot(sector, kPadOuter2SupS);}
  const TGeoHMatrix &PadInner2Tpc(Int_t sector = 1)   const {return TpcRot(sector, kPadInner2Tpc);}
  const TGeoHMatrix &PadOuter2Tpc(Int_t sector = 1)   const {return TpcRot(sector, kPadOuter2Tpc);}
  const TGeoHMatrix &PadInner2Glob(Int_t sector = 1)  const {return TpcRot(sector, kPadInner2Glob);}
  const TGeoHMatrix &PadOuter2Glob(Int_t sector = 1)  const {return TpcRot(sector, kPadOuter2Glob);}

  const TGeoHMatrix &SubS2SupS(Int_t sector = 1, Int_t row = 1) const {Int_t k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kSubSInner2SupS : kSubSOuter2SupS; return TpcRot(sector, k);}
  const TGeoHMatrix &SubS2Tpc(Int_t sector = 1, Int_t row = 1)  const {Int_t k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kSubSInner2Tpc : kSubSOuter2Tpc; return TpcRot(sector, k);}
  const TGeoHMatrix &SubS2Glob(Int_t sector = 1, Int_t row = 1) const {Int_t k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kSubSInner2Glob : kSubSOuter2Glob; return TpcRot(sector, k);}

  const TGeoHMatrix &Pad2SupS(Int_t sector = 1, Int_t row = 1)  const {Int_t k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kPadInner2SupS : kPadOuter2SupS; return TpcRot(sector, k);}
  const TGeoHMatrix &Pad2Tpc(Int_t sector = 1, Int_t row = 1)   const {Int_t k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kPadInner2Tpc : kPadOuter2Tpc; return TpcRot(sector, k);}
  const TGeoHMatrix &Pad2Glob(Int_t sector = 1, Int_t row = 1)  const {Int_t k = (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) ? kPadInner2Glob : kPadOuter2Glob; return TpcRot(sector, k);}
};
#endif
