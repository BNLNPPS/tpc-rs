#pragma once

class Configurator;


class MagField
{
 public:

  enum class MagFieldType {kConstant, kMapped};

  MagField(const tpcrs::Configurator& cfg, MagFieldType field_type = MagFieldType::kMapped, double scale = 1);

  void GetFieldValue(const double x[3], float B[3]) const;
  void GetFieldValue3D(const float x[3], float B[3]) const;

  static void Search(int N, const float Xarray[], float x, int &low);
  float Interpolate(const float Xarray[], const float Yarray[], int order, const float x) const;

 private:

  void ReadField();
  void InterpolateField2D(float r, float z, float &Br_value, float &Bz_value) const;
  void InterpolateField3D(float r, float z, float phi, float &Br_value, float &Bz_value, float &Bphi_value) const;

  enum ESmFSizes {nZ = 57, nR = 28, nPhi = 37};

  const tpcrs::Configurator& cfg_;

  MagFieldType field_type_;

  /// Defined by the value of MagFactor.ScaleFactor in Configurator
  double scale_factor_;

  float Bz[nZ][nR], Br[nZ][nR] ;
  float Radius[nR], ZList[nZ] ;
  float Bz3D[nPhi][nZ][nR], Br3D[nPhi][nZ][nR], Bphi3D[nPhi][nZ][nR] ;
  float R3D[nR], Z3D[nZ], Phi3D[nPhi] ;
};
