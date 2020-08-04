#pragma once

#include "TGeoMatrix.h"
#include "tpcrs/configurator.h"


class MagField
{
 public:

  enum  EBField {kUndefined = 0, kConstant = 1, kMapped = 2};

  MagField (const tpcrs::Configurator& cfg, EBField map = kMapped, float Rescale = 1);

  void BField  (const double x[], float B[]);
  void B3DField(const float x[], float B[]);

  static void Search(int N, const float Xarray[], float x, int &low);
  float Interpolate(const float Xarray[], const float Yarray[], const int ORDER, const float x);

 private:

  void ReadField();
  void Interpolate2DBfield(const float r, const float z, float &Br_value, float &Bz_value);
  void Interpolate3DBfield(const float r, const float z, const float phi, float &Br_value, float &Bz_value, float &Bphi_value);

  enum  ESmFSizes {nZ = 57, nR = 28, nPhi = 37};

  const tpcrs::Configurator& cfg_;

  TGeoRotation fMagFieldRotation;

  EBField  fMap;     // (D) = kMapped; Global flag to indicate static arrays are full
  float  fFactor;    // (D) = 1.0    ; Multiplicative factor (allows scaling and sign reversal)
  float  fRescale;   // (D) = 1.0    ; Multiplicative factor (allows re-scaling wrt which map read)

  float  Bz[nZ][nR], Br[nZ][nR] ;
  float  Radius[nR], ZList[nZ] ;
  float  Bz3D[nPhi][nZ][nR], Br3D[nPhi][nZ][nR], Bphi3D[nPhi][nZ][nR] ;
  float  R3D[nR], Z3D[nZ], Phi3D[nPhi] ;
};
