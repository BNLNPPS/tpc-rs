/***************************************************************************
 * Author:  David Hardtke
 *
 * Description: This is the interface between the database and the offline
 *              TPC software.  This classes takes care of the annoying
 *              calls to the root infrastucture, packages and manipulates
 *              the data, and returns the data to the user via simple
 *              interface classes.
 **************************************************************************/

#include "TEnv.h"
#include "TGeoManager.h"
#include "TVector3.h"
#include "TString.h"

#include "enums.h"
#include "logger.h"
#include "mag_field.h"
#include "tpc_db.h"

bool StTpcDb::mOldScheme = true;


StTpcDb::StTpcDb()
{
  mTpc2GlobMatrix = new TGeoHMatrix("Default Tpc2Glob");

  for (int i = 1; i <= 24; i++) {
    for (int k = 0; k < kTotalTpcSectorRotaions; k++) {
      mTpcSectorRotations[i - 1][k] = new TGeoHMatrix(Form("Default %02i %i", i, k));
    }
  }

  mFlip = new TGeoHMatrix;
  mzGG = St_tpcDimensionsC::instance()->gatingGridZ(); // zGG
  double Rotation[9] = {0, 1, 0,
                          1, 0, 0,
                          0, 0, -1
                         };
  //  double Translation[3] = {0, 0, mzGG};
  mFlip->SetName("Flip"); mFlip->SetRotation(Rotation);// mFlip->SetTranslation(Translation);
  mSwap[0] = new TGeoTranslation("Signed Drift distance to z for East", 0, 0, -mzGG);
  mSwap[1] = new TGeoTranslation("Signed Drift distance to z for West", 0, 0,  mzGG);
  mHalf[0] = new TGeoHMatrix("Default for east part of TPC");
  mHalf[1] = new TGeoHMatrix("Default for west part of TPC");

  SetDriftVelocity();
  SetTpcRotations();
}




StTpcDb::~StTpcDb()
{
  for (int i = 0; i < 24; i++) {
    for (int k = 0; k < kTotalTpcSectorRotaions; k++)
      SafeDelete(mTpcSectorRotations[i][k]);
  }

  SafeDelete(mHalf[0]);
  SafeDelete(mHalf[1]);
  SafeDelete(mSwap[0]);
  SafeDelete(mSwap[1]);
  SafeDelete(mTpc2GlobMatrix);
  SafeDelete(mFlip);
}

float StTpcDb::DriftVelocity(int sector)
{
  TPC::Half half = (sector <= 12 ? TPC::Half::first : TPC::Half::second);
  return 1e6 * mDriftVel[half];
}

void StTpcDb::SetDriftVelocity()
{
  tpcDriftVelocity_st* dv = St_tpcDriftVelocityC::instance()->Struct();
  mDriftVel[0] = dv->laserDriftVelocityWest > 0 ? dv->laserDriftVelocityWest : dv->cathodeDriftVelocityWest;
  mDriftVel[1] = dv->laserDriftVelocityEast > 0 ? dv->laserDriftVelocityEast : dv->cathodeDriftVelocityEast;
}


void StTpcDb::SetTpcRotations()
{
  // Pad [== sector12 == localsector (SecL, ideal)] => subsector (SubS,local sector aligned) => flip => sector (SupS) => tpc => global
  //                                    ------
  //old:  global = Tpc2GlobalMatrix() * SupS2Tpc(sector) *                                    Flip() * {SubSInner2SupS(sector) | SubSOuter2SupS(sector)}
  //new:  global = Tpc2GlobalMatrix() * SupS2Tpc(sector) * StTpcSuperSectorPosition(sector) * Flip() * {                     I | StTpcOuterSectorPosition(sector)}
  //      StTpcSuperSectorPosition(sector) * Flip() = Flip() * SubSInner2SupS(sector)
  // =>  StTpcSuperSectorPosition(sector) = Flip() * SubSInner2SupS(sector) * Flip()^-1
  //      StTpcSuperSectorPosition(sector) * Flip() * StTpcOuterSectorPosition(sector) = Flip() *  SubSOuter2SupS(sector)
  // =>  StTpcOuterSectorPosition(sector) = Flip()^-1 * StTpcSuperSectorPosition(sector)^-1 *  Flip() *  SubSOuter2SupS(sector)
  /*
    .                                                                                             <-- the system of coordinate where Outer to Inner Alignment done -->
    global = Tpc2GlobalMatrix() * SupS2Tpc(sector) * StTpcSuperSectorPosition(sector) * Flip() * {                     I | StTpcOuterSectorPosition(sector)} * local
    .                                                result of super sector alignment                                      result of Outer to Inner sub sector alignment
  */
  /* 03/07/14
     StTpcPadCoordinate P(sector,row,pad,timebacket);
     x = xFromPad()
     z = zFromTB() - drift distance
     StTpcLocalSectorCoordinate LS(position(x,y,z),sector,row);

     StTpcLocalCoordinate  Tpc(position(x,y,z),sector,row);
     Pad2Tpc(sector,row).LocalToMaster(LS.position().xyz(), Tpc.postion().xyz())
     Flip transformation from Pad Coordinate system (xFromPad(pad), yFromRow(row), DriftDistance(timebacket)) => (y, x, -z): local sector CS => super sectoe CS



             (0 1  0 0  ) ( x )    ( y )
     Flip ;  (1 0  0 0  ) ( y ) =  ( x )
             (0 0 -1 zGG) ( z )    ( zGG - z)
             (0 0  0 1  ) ( 1 )    ( 1 )
     Z_tpc is not changed during any sector transformation  !!!


   */
  //  TGeoTranslation T123(0,123,0); T123.SetName("T123"); if (Debug() > 1) T123.Print();
  assert(St_tpcDimensionsC::instance()->numberOfSectors() == 24);
  float gFactor = St_MagFactorC::instance()->ScaleFactor();
  double phi, theta, psi;
  int iphi;
  TGeoRotation* rotm = 0;
  TObjArray* listOfMatrices = 0;
  TString Rot;

  if (gEnv->GetValue("NewTpcAlignment", 0) != 0) mOldScheme = false;

  if (! mOldScheme) {
    LOG_INFO << "StTpcDb::SetTpcRotations use new schema for Rotation matrices\n";
  }
  else {
    LOG_INFO << "StTpcDb::SetTpcRotations use old schema for Rotation matrices\n";
  }

  for (int sector = 0; sector <= 24; sector++) {// loop over Tpc as whole, sectors, inner and outer subsectors
    int k;
    int k1 = kSupS2Tpc;
    int k2 = kTotalTpcSectorRotaions;

    if (sector == 0) {k2 = k1; k1 = kUndefSector;}

    for (k = k1; k < k2; k++) {
      int Id     = 0;
      TGeoHMatrix rotA; // After alignment

      if (!sector ) { // TPC Reference System
        if (mOldScheme) { // old scheme
          St_tpcGlobalPositionC* tpcGlobalPosition = St_tpcGlobalPositionC::instance();
          assert(tpcGlobalPosition);
          Id = 1;
          phi   = 0.0;                                               // -gamma large uncertainty, so set to 0
          theta = tpcGlobalPosition->PhiXZ_geom() * 180.0 / M_PI; // -beta
          psi   = tpcGlobalPosition->PhiYZ_geom() * 180.0 / M_PI; // -alpha
          rotA.RotateX(-psi);
          rotA.RotateY(-theta);
          rotA.RotateZ(-phi);
          double transTpcRefSys[3] = {tpcGlobalPosition->LocalxShift(),
                                        tpcGlobalPosition->LocalyShift(),
                                        tpcGlobalPosition->LocalzShift()
                                       };
          rotA.SetTranslation(transTpcRefSys);
        }
        else {
          rotA = StTpcPosition::instance()->GetMatrix();
          *mHalf[TPC::Half::first] = StTpcHalfPosition::instance()->GetEastMatrix();
          *mHalf[TPC::Half::second] = StTpcHalfPosition::instance()->GetWestMatrix();
        }
      }
      else {
        Id = 10 * sector + k;
        TPC::Half part = TPC::Half::first;

        if (sector <= 12) part = TPC::Half::second;

        switch (k) {
        case kSupS2Tpc: // SupS => Tpc
          if (sector <= 12) {iphi = (360 + 90 - 30 * sector      ) % 360; Rot = Form("R%03i", iphi);}
          else              {iphi = (      90 + 30 * (sector - 12)) % 360; Rot = Form("Y%03i", iphi);}

          rotm = 0;

          if (gGeoManager) {
            listOfMatrices =  gGeoManager->GetListOfMatrices();
            rotm = (TGeoRotation*) listOfMatrices->FindObject(Rot);
          }

          if (! rotm) {
            if (sector <= 12) rotm = new TGeoRotation(Rot);
            else              rotm = new TGeoRotation(Rot,   90.0,    0.0,  90.0,  -90.0,  180.0,    0.00); // Flip (x,y,z) => ( x,-y,-z)

            rotm->RotateZ(iphi);
          }

          rotA = (*mSwap[part]) * (*mHalf[part]) * (*rotm);
          rotA *= StTpcSuperSectorPosition::instance()->GetMatrix(sector - 1);

          if (gGeoManager) rotm->RegisterYourself();
          else             SafeDelete(rotm);

          break;

        case kSupS2Glob:      // SupS => Tpc => Glob
          rotA = Tpc2GlobalMatrix() * SupS2Tpc(sector);
          break;

        case kSubSInner2SupS:
          if (mOldScheme) 	  rotA = Flip();
          else                    rotA = Flip() * StTpcInnerSectorPosition::instance()->GetMatrix(sector - 1);

          break;

        case kSubSOuter2SupS:
          if (mOldScheme) rotA = Flip() * StTpcOuterSectorPosition::instance()->GetMatrix(sector - 1);
          else           {
            rotA = Flip() * StTpcOuterSectorPosition::instance()->GetMatrix(sector - 1);

            if (StTpcOuterSectorPosition::instance()->GetNRows() > 24) {
              if (gFactor > 0.2) {
                rotA *= StTpcOuterSectorPosition::instance()->GetMatrix(sector - 1 + 24);
              }
              else if (gFactor < -0.2) {
                rotA *= StTpcOuterSectorPosition::instance()->GetMatrix(sector - 1 + 24).Inverse();
              }
            }
          }

          break;

        case kSubSInner2Tpc:  rotA = SupS2Tpc(sector) * SubSInner2SupS(sector); break; // (Subs[io] => SupS) => Tpc

        case kSubSOuter2Tpc:  rotA = SupS2Tpc(sector) * SubSOuter2SupS(sector); break; // -"-

        case kSubSInner2Glob: rotA = Tpc2GlobalMatrix() * SubSInner2Tpc(sector);  break; // Subs[io] => SupS => Tpc) => Glob

        case kSubSOuter2Glob: rotA = Tpc2GlobalMatrix() * SubSOuter2Tpc(sector);  break; // -"-

        case kPadInner2SupS:  rotA = SubSInner2SupS(sector); break; // (Pad == SecL) => (SubS[io] => SupS)

        case kPadOuter2SupS:  rotA = SubSOuter2SupS(sector); break; // -"-

        case kPadInner2Tpc:   rotA = SupS2Tpc(sector) * PadInner2SupS(sector); break; // (Pad == SecL) => (SubS[io] => SupS => Tpc)

        case kPadOuter2Tpc:   rotA = SupS2Tpc(sector) * PadOuter2SupS(sector); break; // -"-

        case kPadInner2Glob:  rotA = Tpc2GlobalMatrix() * PadInner2Tpc(sector); break; // (Pad == SecL) => (SubS[io] => SupS => Tpc => Glob)

        case kPadOuter2Glob:  rotA = Tpc2GlobalMatrix() * PadOuter2Tpc(sector); break; // -"-

        default:
          assert(0);
        }
      }

      // Normalize
      double* r = rotA.GetRotationMatrix();
      double norm;
      TVector3 d(r[0], r[3], r[6]); norm = 1 / d.Mag(); d *= norm;
      TVector3 t(r[2], r[5], r[8]); norm = 1 / t.Mag(); t *= norm;
      TVector3 n(r[1], r[4], r[7]);
      TVector3 c = d.Cross(t);

      if (c.Dot(n) < 0) c *= -1;

      double rot[9] = {
        d[0], c[0], t[0],
        d[1], c[1], t[1],
        d[2], c[2], t[2]
      };
      rotA.SetRotation(rot);
      const char* names[kTotalTpcSectorRotaions] = {
        "SupS_%02itoTpc",
        "SupS_%02itoGlob",
        "SubS_%02iInner2SupS",
        "SubS_%02iOuter2SupS",
        "SubS_%02iInner2Tpc",
        "SubS_%02iOuter2Tpc",
        "SubS_%02iInner2Glob",
        "SubS_%02iOuter2Glob",
        "PadInner2SupS_%02i",
        "PadOuter2SupS_%02i",
        "SupS_%02i12Inner2Tpc",
        "SupS_%02i12Outer2Tpc",
        "SupS_%02i12Inner2Glob",
        "SupS_%02i12Outer2Glob"
      };

      if (sector == 0) rotA.SetName("Tpc2Glob");
      else             rotA.SetName(Form(names[k], sector));

      SetTpcRotationMatrix(&rotA, sector, k);
    }
  }
}
