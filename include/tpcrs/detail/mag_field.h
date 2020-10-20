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
    if (z >= c2_.front() && z <= c2_.back() && r <= c1_.back())
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

  static float Interpolate(const float xs[2], const float ys[2], const float x);

 private:

  void ReadValues();
  void InterpolateField2D(double r, double z, double &Br_value, double &Bz_value) const;

  const tpcrs::Configurator& cfg_;

  MagFieldType field_type_;

  /// Defined by the value of MagFactor.ScaleFactor in Configurator
  double scale_factor_;

  std::vector<float> c1_;
  std::vector<float> c2_;
  std::vector<float> c3_;

  std::vector<float> v1_;
  std::vector<float> v2_;
  std::vector<float> v3_;
};

} }
