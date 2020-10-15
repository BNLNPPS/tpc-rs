#pragma once


namespace tpcrs { namespace detail {

class Configurator;


class MagField
{
 public:

  enum class MagFieldType {kConstant, kMapped};

  MagField(const tpcrs::Configurator& cfg, MagFieldType field_type, double scale);

  MagField(const tpcrs::Configurator& cfg) : MagField(cfg, MagFieldType::kMapped, cfg.S<MagFactor>().ScaleFactor) {}


  /**
   * Returns the value of magnetic field in Cartesian coordinates
   */
  template<typename Vec3>
  Vec3 ValueAt(Vec3 p) const
  {
    Vec3 B{};

    double B_r = 0;
    double B_z = 0;
    double r = std::sqrt(p.x*p.x + p.y*p.y);
    double z = p.z;

    // within map
    if (z >= c2_[0] && z <= c2_[nZ - 1] && r <= c1_[nR - 1])
    {
      InterpolateField2D(r, z, B_r, B_z);

      if (r != 0) {
        B.x = B_r * (p.x / r) ;
        B.y = B_r * (p.y / r) ;
      }

      B.z = B_z;
    }

    return B;
  }

  static void Search(const int N, const float Xarray[], const float x, int &low);
  static float Interpolate(const float xs[2], const float ys[2], const float x);

 private:

  void ReadValues();
  void InterpolateField2D(double r, double z, double &Br_value, double &Bz_value) const;
  void InterpolateField3D(float r, float z, float phi, float &Br_value, float &Bz_value, float &Bphi_value) const;

  enum ESmFSizes {nR = 28, nZ = 57, nPhi = 37};

  const tpcrs::Configurator& cfg_;

  MagFieldType field_type_;

  /// Defined by the value of MagFactor.ScaleFactor in Configurator
  double scale_factor_;

  float R_[nR], Z_[nZ], Bz[nZ][nR], Br[nZ][nR];

  float R3D[nR], Z3D[nZ], Phi3D[nPhi], Bz3D[nPhi][nZ][nR], Br3D[nPhi][nZ][nR], Bphi3D[nPhi][nZ][nR];


  std::vector<float> c1_;
  std::vector<float> c2_;
  std::vector<float> c3_;

  std::vector<float> v1_;
  std::vector<float> v2_;
  std::vector<float> v3_;
};

} }
