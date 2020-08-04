#include "TEnv.h"
#include "TGeoManager.h"
#include "TVector3.h"
#include "TString.h"

#include "tpcrs/configurator.h"
#include "coords.h"
#include "logger.h"
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


CoordTransform::CoordTransform(const tpcrs::Configurator& cfg) :
  cfg_(cfg)
{
  mTimeBinWidth = 1e6 / cfg_.S<starClockOnl>().frequency;
  mInnerSectorzOffset = cfg_.S<tpcEffectiveGeom>().z_inner_offset;
  mOuterSectorzOffset = cfg_.S<tpcEffectiveGeom>().z_outer_offset;

  mTpc2GlobMatrix = new TGeoHMatrix("Default Tpc2Glob");

  for (int i = 1; i <= 24; i++) {
    for (int k = 0; k < kTotalTpcSectorRotaions; k++) {
      mTpcSectorRotations[i - 1][k] = new TGeoHMatrix(Form("Default %02i %i", i, k));
    }
  }

  mFlip = new TGeoHMatrix;
  mzGG = cfg_.S<tpcPadPlanes>().outerSectorPadPlaneZ - cfg_.S<tpcWirePlanes>().outerSectorGatingGridPadSep;

  double Rotation[9] = {0, 1, 0,
                        1, 0, 0,
                        0, 0, -1};

  mFlip->SetName("Flip"); mFlip->SetRotation(Rotation);// mFlip->SetTranslation(Translation);
  mSwap[0] = new TGeoTranslation("Signed Drift distance to z for East", 0, 0, -mzGG);
  mSwap[1] = new TGeoTranslation("Signed Drift distance to z for West", 0, 0,  mzGG);
  mHalf[0] = new TGeoHMatrix("Default for east part of TPC");
  mHalf[1] = new TGeoHMatrix("Default for west part of TPC");

  SetTpcRotations();
}


// Local Sector Coordnate <-> Tpc Raw Pad Coordinate
void CoordTransform::local_sector_to_hardware(const StTpcLocalSectorCoordinate &a, StTpcPadCoordinate &b, bool useT0, bool useTau) const
{
  // useT0 = true for pad and false for cluster, useTau = true for data cluster and  = false for MC
  int row = a.row;

  if (row < 1 || row > cfg_.S<tpcPadPlanes>().padRows)
    row = rowFromLocalY(a.position.y, a.sector);

  double probablePad = padFromX(a.position.x, a.sector, a.row);
  double zoffset = (row > cfg_.S<tpcPadPlanes>().innerPadRows) ? mOuterSectorzOffset : mInnerSectorzOffset;
  double t0offset = (useT0 && a.sector >= 1 && a.sector <= 24) ? cfg_.S<tpcPadGainT0>().T0[a.sector-1][row-1][tpcrs::irint(probablePad)-1] : 0;
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
  double zoffset =  (a.row > cfg_.S<tpcPadPlanes>().innerPadRows) ? mOuterSectorzOffset : mInnerSectorzOffset;
  double t0offset = useT0 ? cfg_.S<tpcPadGainT0>().T0[a.sector-1][a.row-1][tpcrs::irint(a.pad)-1] : 0;
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
  if (row > cfg_.S<tpcPadPlanes>().padRows) row = cfg_.S<tpcPadPlanes>().padRows;

  double pitch = (row <= cfg_.S<tpcPadPlanes>().innerPadRows) ?
                         cfg_.S<tpcPadPlanes>().innerSectorPadPitch :
                         cfg_.S<tpcPadPlanes>().outerSectorPadPitch;

  int npads = cfg_.C<St_tpcPadPlanesC>().numberOfPadsAtRow(row);
  double probablePad = (npads + 1.) / 2. - x / pitch;

  // CAUTION: pad cannot be <1
  if (probablePad < 0.500001) {
    probablePad = 0.500001;
  }

  return probablePad;
}


double CoordTransform::xFromPad(int sector, int row, double pad) const      // x coordinate in sector 12
{
  if (row > cfg_.S<tpcPadPlanes>().padRows) row = cfg_.S<tpcPadPlanes>().padRows;

  double pitch = (row <= cfg_.S<tpcPadPlanes>().innerPadRows) ?
                         cfg_.S<tpcPadPlanes>().innerSectorPadPitch :
                         cfg_.S<tpcPadPlanes>().outerSectorPadPitch;

  int npads = cfg_.C<St_tpcPadPlanesC>().numberOfPadsAtRow(row);
  double xPad = -pitch * (pad - (npads + 1.) / 2.);

  return xPad;
}
// Coordinate from Row
//
//Local Transformation...


double CoordTransform::zFromTB(double tb, int sector, int row, int pad) const
{
  if (row > cfg_.S<tpcPadPlanes>().padRows) row = cfg_.S<tpcPadPlanes>().padRows;

  double trigT0 = triggerTimeOffset() * 1e6; // units are s
  double elecT0 = cfg_.S<tpcElectronics>().tZero;    // units are us
  double sectT0 = cfg_.S<tpcPadrowT0>(sector-1).T0[row-1];  // units are us
  double t0 = trigT0 + elecT0 + sectT0;
  int l = sector;

  if (tpcrs::IsInner(row, cfg_)) l += 24;

  double tbx = tb + cfg_.S<tpcSectorT0offset>().t0[l-1];

  double time = t0 + tbx * mTimeBinWidth;
  double z = tpcrs::DriftVelocity(sector, cfg_) * 1e-6 * time;
  return z;
}


double CoordTransform::tBFromZ(double z, int sector, int row, int pad) const
{
  if (row > cfg_.S<tpcPadPlanes>().padRows) row = cfg_.S<tpcPadPlanes>().padRows;

  double trigT0 = triggerTimeOffset() * 1e6; // units are s
  double elecT0 = cfg_.S<tpcElectronics>().tZero;    // units are us
  double sectT0 = cfg_.S<tpcPadrowT0>(sector-1).T0[row-1];  // units are us
  double t0 = trigT0 + elecT0 + sectT0;
  double time = z / (tpcrs::DriftVelocity(sector, cfg_) * 1e-6);
  int l = sector;

  if (tpcrs::IsInner(row, cfg_)) l += 24;

  double tb = (time - t0) / mTimeBinWidth - cfg_.S<tpcSectorT0offset>().t0[l-1];

  return tb;
}


int CoordTransform::rowFromLocalY(double y, int sector) const
{
  static int Nrows = 0;
  static double* Radii = 0;

  if (Nrows != cfg_.S<tpcPadPlanes>().padRows) {
    Nrows = cfg_.S<tpcPadPlanes>().padRows;

    if (Radii) delete [] Radii;

    Radii = new double[Nrows + 1];

    for (int i = 1; i <= Nrows + 1; i++) {
      if (i == 1) {
        Radii[i - 1] =  (3 * tpcrs::RadialDistanceAtRow(i, cfg_)
                           - tpcrs::RadialDistanceAtRow(i + 1, cfg_)) / 2;
      }
      else if (i == Nrows + 1) {
        Radii[i - 1] =  (3 * tpcrs::RadialDistanceAtRow(i - 1, cfg_)
                           - tpcrs::RadialDistanceAtRow(i - 2, cfg_)) / 2;
      }
      else {
        Radii[i - 1] = (tpcrs::RadialDistanceAtRow(i - 1, cfg_) +
                        tpcrs::RadialDistanceAtRow(i, cfg_)) / 2;
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

  if (row < 1 || row > cfg_.S<tpcPadPlanes>().padRows)
    row = rowFromLocalY(a.position.y, a.sector);

  Coords xGG;
  Pad2Tpc(a.sector, row).LocalToMasterVect(a.position.xyz(), xGG.xyz());
  const double* trans = Pad2Tpc(a.sector, row).GetTranslation();
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  GG2TPC.LocalToMaster(xGG.xyz(), b.position.xyz());

  b.row = row;
  b.sector = a.sector;
}


void CoordTransform::local_to_local_sector(const StTpcLocalCoordinate &a, StTpcLocalSectorCoordinate &b) const
{
  int row = a.row;

  if (row < 1 || row > cfg_.S<tpcPadPlanes>().padRows) {
    Coords xyzS;
    SupS2Tpc(a.sector).MasterToLocalVect(a.position.xyz(), xyzS.xyz());
    row = rowFromLocalY(xyzS.x, a.sector);
  }

  const double* trans = Pad2Tpc(a.sector, row).GetTranslation();
  TGeoTranslation GG2TPC(trans[0], trans[1], trans[2]);
  Coords xGG;
  GG2TPC.MasterToLocal(a.position.xyz(), xGG.xyz());
  Pad2Tpc(a.sector, row).MasterToLocalVect(xGG.xyz(), b.position.xyz());

  b.row = row;
  b.sector = a.sector;
}

bool CoordTransform::mOldScheme = true;

CoordTransform::~CoordTransform()
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


void CoordTransform::SetTpcRotations()
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
  assert(cfg_.S<tpcDimensions>().numberOfSectors == 24);
  float gFactor = cfg_.S<MagFactor>().ScaleFactor;
  double phi, theta, psi;
  int iphi;
  TGeoRotation* rotm = 0;
  TObjArray* listOfMatrices = 0;
  TString Rot;

  if (gEnv->GetValue("NewTpcAlignment", 0) != 0) mOldScheme = false;

  for (int sector = 0; sector <= 24; sector++) {// loop over Tpc as whole, sectors, inner and outer subsectors
    int k;
    int k1 = kSupS2Tpc;
    int k2 = kTotalTpcSectorRotaions;

    if (sector == 0) {k2 = k1; k1 = kUndefSector;}

    for (k = k1; k < k2; k++) {
      int Id = 0;
      TGeoHMatrix rotA; // After alignment

      if (!sector) { // TPC Reference System
        if (mOldScheme) { // old scheme
          St_tpcGlobalPositionC& tpcGlobalPosition = cfg_.C<St_tpcGlobalPositionC>();
          Id = 1;
          phi   = 0.0;                                               // -gamma large uncertainty, so set to 0
          theta = tpcGlobalPosition.PhiXZ_geom() * 180.0 / M_PI; // -beta
          psi   = tpcGlobalPosition.PhiYZ_geom() * 180.0 / M_PI; // -alpha
          rotA.RotateX(-psi);
          rotA.RotateY(-theta);
          rotA.RotateZ(-phi);
          double transTpcRefSys[3] = {tpcGlobalPosition.LocalxShift(),
                                      tpcGlobalPosition.LocalyShift(),
                                      tpcGlobalPosition.LocalzShift() };
          rotA.SetTranslation(transTpcRefSys);
        }
        else {
          rotA = cfg_.C<StTpcPosition>().GetMatrix();
          *mHalf[tpcrs::TPC::Half::first]  = cfg_.C<StTpcHalfPosition>().GetEastMatrix();
          *mHalf[tpcrs::TPC::Half::second] = cfg_.C<StTpcHalfPosition>().GetWestMatrix();
        }
      }
      else {
        Id = 10 * sector + k;
        tpcrs::TPC::Half part = tpcrs::TPC::Half::first;

        if (sector <= 12) part = tpcrs::TPC::Half::second;

        switch (k) {
        case kSupS2Tpc: // SupS => Tpc
          if (sector <= 12) {iphi = (360 + 90 - 30 * sector      ) % 360; Rot = Form("R%03i", iphi);}
          else              {iphi = (      90 + 30 * (sector - 12)) % 360; Rot = Form("Y%03i", iphi);}

          rotm = 0;

          if (gGeoManager) {
            listOfMatrices =  gGeoManager->GetListOfMatrices();
            rotm = (TGeoRotation*) listOfMatrices->FindObject(Rot);
          }

          if (!rotm) {
            if (sector <= 12) rotm = new TGeoRotation(Rot);
            else              rotm = new TGeoRotation(Rot,   90.0,    0.0,  90.0,  -90.0,  180.0,    0.00); // Flip (x,y,z) => ( x,-y,-z)

            rotm->RotateZ(iphi);
          }

          rotA = (*mSwap[part]) * (*mHalf[part]) * (*rotm);
          rotA *= cfg_.C<StTpcSuperSectorPosition>().GetMatrix(sector - 1);

          if (gGeoManager) rotm->RegisterYourself();
          else             SafeDelete(rotm);

          break;

        case kSupS2Glob:      // SupS => Tpc => Glob
          rotA = Tpc2GlobalMatrix() * SupS2Tpc(sector);
          break;

        case kSubSInner2SupS:
          if (mOldScheme) rotA = Flip();
          else            rotA = Flip() * cfg_.C<StTpcInnerSectorPosition>().GetMatrix(sector - 1);

          break;

        case kSubSOuter2SupS:
          if (mOldScheme)
            rotA = Flip() * cfg_.C<StTpcOuterSectorPosition>().GetMatrix(sector - 1);
          else
          {
            rotA = Flip() * cfg_.C<StTpcOuterSectorPosition>().GetMatrix(sector - 1);

            if (cfg_.C<StTpcOuterSectorPosition>().GetNRows() > 24) {
              if (gFactor > 0.2) {
                rotA *= cfg_.C<StTpcOuterSectorPosition>().GetMatrix(sector - 1 + 24);
              }
              else if (gFactor < -0.2) {
                rotA *= cfg_.C<StTpcOuterSectorPosition>().GetMatrix(sector - 1 + 24).Inverse();
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
