/***********************************************************************
 * Author: Jim Thomas   11/1/2000
 *
 * Description: the STAR Magnetic Field
  ***********************************************************************/

#ifndef TPCRS_MAG_FIELD_H_
#define TPCRS_MAG_FIELD_H_

#include <cmath>

#include "TH2.h"
#include "TGeoMatrix.h"


class StarMagField
{
 public:
  enum   EBField  { kUndefined = 0, kConstant = 1, kMapped = 2, kChain = 3 } ;
  enum   ESmFSizes {nZ = 57, nR = 28, nPhi = 37, nZSteel = 16, nRSteel = 115, nPhiSteel = 25};
  static  void    Search ( int N, const float Xarray[], float x, int &low ) ;
  float Interpolate ( const float Xarray[], const float Yarray[],
                                const int ORDER, const float x ) ;
  void    Interpolate2DBfield ( const float r, const float z,
                                        float &Br_value, float &Bz_value ) ;
  void    Interpolate2ExtDBfield ( const float r, const float z,
      float &Br_value, float &Bz_value ) ;
  void    Interpolate3DBfield ( const float r, const float z, const float phi,
                                        float &Br_value, float &Bz_value, float &Bphi_value ) ;
 private:

  void    ReadField ( ) ;
  TGeoRotation fStarMagFieldRotation;
  TH2F* fBzdZCorrection; // correction due to endcap calomiter
  TH2F* fBrdZCorrection; // correction due to endcap calomiter

 public:

  void    Interpolate3DBSteelfield ( const float r, const float z, const float phi,
      float &Br_value, float &Bz_value, float &Bphi_value );

  EBField  fMap;       // (D) = kMapped; Global flag to indicate static arrays are full
  float  fFactor;    // (D) = 1.0    ; Multiplicative factor (allows scaling and sign reversal)
  float  fRescale;   // (D) = 1.0    ; Multiplicative factor (allows re-scaling wrt which map read)
  float  fBDipole;   // (D) = -42.67 ; field value (kG)
  float  fRmaxDip;   // (D) =  15.34 ; Inner field volume radius
  float  fZminDip;   // (D) =  980.0 ; StArt of the DX mAgnet in Z
  float  fZmaxDip;   // (D) = 1350.0 ; End of the DX mAgnet in Z
  bool   fLock;      // (D) = false ; Set true if lock above values

  float  Bz[nZ][nR], Br[nZ][nR] ;
  float  Radius[nR], ZList[nZ] ;
  float  Bz3D[nPhi][nZ][nR], Br3D[nPhi][nZ][nR], Bphi3D[nPhi][nZ][nR] ;
  float  R3D[nR], Z3D[nZ], Phi3D[nPhi] ;
  float  R3DSteel[nRSteel], Z3DSteel[nZSteel], Phi3DSteel[nPhiSteel] ;
  float  Bz3DSteel[nPhiSteel][nZSteel][nRSteel];
  float  Bx3DSteel[nPhiSteel][nZSteel][nRSteel], By3DSteel[nPhiSteel][nZSteel][nRSteel] ;
  static bool     mConstBz;

  float  Br3DSteel[nPhiSteel][nZSteel][nRSteel], Bphi3DSteel[nPhiSteel][nZSteel][nRSteel] ;

 private:

  StarMagField ( EBField map     = kMapped, float Factor  =      1,
                 bool  Lock    =  false, float Rescale =      1,
                 float Bdipole =  -42.67, float Rmaxdip =  15.34,
                 float Zmindip =   980.0, float Zmaxdip = 1350.0) ;
  ~StarMagField ()
  {
    SafeDelete(fBzdZCorrection);
    SafeDelete(fBrdZCorrection);
  }
 public:

  static StarMagField& Instance()
  {
    static StarMagField instance;
    return instance;
  }

  void    BField   ( const float x[], float B[] ) ;
  void    BField   ( const double x[], double B[] ) ;
  void    B3DField ( const float x[], float B[] ) ;
  const TGeoRotation &StarMagFieldRotation() {return * &fStarMagFieldRotation;}
};

#endif
