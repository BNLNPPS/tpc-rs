#pragma once

#include "TH2.h"
#include "TGeoMatrix.h"
#include "tpcrs/configurator.h"


class MagField
{
 public:

  enum  EBField {kUndefined = 0, kConstant = 1, kMapped = 2};

  MagField (const tpcrs::Configurator& cfg, EBField map     = kMapped,
                 float Rescale =      1,
                 float Bdipole =  -42.67, float Rmaxdip =  15.34,
                 float Zmindip =   980.0, float Zmaxdip = 1350.0) ;

  ~MagField ()
  {
    SafeDelete(fBzdZCorrection);
    SafeDelete(fBrdZCorrection);
  }

  void BField  (const double x[], float B[]);
  void B3DField(const float x[], float B[]);

  static void Search(int N, const float Xarray[], float x, int &low);
  float Interpolate(const float Xarray[], const float Yarray[], const int ORDER, const float x);

 private:

  void ReadField();
  void Interpolate2DBfield(const float r, const float z, float &Br_value, float &Bz_value);
  void  Interpolate2ExtDBfield ( const float r, const float z, float &Br_value, float &Bz_value ) ;
  void Interpolate3DBfield(const float r, const float z, const float phi, float &Br_value, float &Bz_value, float &Bphi_value);

  enum  ESmFSizes {nZ = 57, nR = 28, nPhi = 37, nZSteel = 16, nRSteel = 115, nPhiSteel = 25};

  const tpcrs::Configurator& cfg_;

  TGeoRotation fMagFieldRotation;
  TH2F* fBzdZCorrection; // correction due to endcap calomiter
  TH2F* fBrdZCorrection; // correction due to endcap calomiter

  void  Interpolate3DBSteelfield ( const float r, const float z, const float phi, float &Br_value, float &Bz_value, float &Bphi_value );

  EBField  fMap;     // (D) = kMapped; Global flag to indicate static arrays are full
  float  fFactor;    // (D) = 1.0    ; Multiplicative factor (allows scaling and sign reversal)
  float  fRescale;   // (D) = 1.0    ; Multiplicative factor (allows re-scaling wrt which map read)
  float  fBDipole;   // (D) = -42.67 ; field value (kG)
  float  fRmaxDip;   // (D) =  15.34 ; Inner field volume radius
  float  fZminDip;   // (D) =  980.0 ; Start of the DX magnet in Z
  float  fZmaxDip;   // (D) = 1350.0 ; End of the DX magnet in Z

  float  Bz[nZ][nR], Br[nZ][nR] ;
  float  Radius[nR], ZList[nZ] ;
  float  Bz3D[nPhi][nZ][nR], Br3D[nPhi][nZ][nR], Bphi3D[nPhi][nZ][nR] ;
  float  R3D[nR], Z3D[nZ], Phi3D[nPhi] ;
  float  R3DSteel[nRSteel], Z3DSteel[nZSteel], Phi3DSteel[nPhiSteel] ;
  float  Bz3DSteel[nPhiSteel][nZSteel][nRSteel];
  float  Bx3DSteel[nPhiSteel][nZSteel][nRSteel], By3DSteel[nPhiSteel][nZSteel][nRSteel] ;

  float  Br3DSteel[nPhiSteel][nZSteel][nRSteel], Bphi3DSteel[nPhiSteel][nZSteel][nRSteel] ;
};
