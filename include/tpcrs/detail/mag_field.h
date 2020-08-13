#pragma once


namespace tpcrs { namespace detail {

class Configurator;


class MagField
{
 public:

  enum class MagFieldType {kConstant, kMapped};

  MagField(const tpcrs::Configurator& cfg, MagFieldType field_type = MagFieldType::kMapped, double scale = 1);


  /**
   * Returns the value of magnetic field in Cartesian coordinates
   */
  template<typename Vec3>
  Vec3 ValueAt(Vec3 p) const
  {
    Vec3 B{};

    double Br_value = 0;
    double Bz_value = 0;
    double r = std::sqrt(p.x*p.x + p.y*p.y);
    double z = p.z;

    // within map
    if (z >= ZList[0] && z <= ZList[nZ - 1] && r <= Radius[nR - 1])
    {
      InterpolateField2D(r, z, Br_value, Bz_value);

      if (r != 0) {
        B.x = Br_value * (p.x / r) ;
        B.y = Br_value * (p.y / r) ;
      }

      B.z = Bz_value;
    }

    return B;
  }

  void GetFieldValue3D(const float x[3], float B[3]) const;

  static void Search(int N, const float Xarray[], float x, int &low);
  float Interpolate(const float Xarray[], const float Yarray[], int order, const float x) const;

 private:

  void ReadField();
  void InterpolateField2D(double r, double z, double &Br_value, double &Bz_value) const;
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

} }
