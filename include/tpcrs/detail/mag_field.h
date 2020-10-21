#pragma once


namespace tpcrs { namespace detail {

class Configurator;


class MagField
{
 public:

  MagField(const tpcrs::Configurator& cfg);

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

  std::vector<float> c1_;
  std::vector<float> c2_;
  std::vector<float> c3_;

  std::vector<float> v1_;
  std::vector<float> v2_;
  std::vector<float> v3_;
};


template<typename T>
struct Interpolator
{
  /**
   * Linear interpolation at point p
   *
   *        (vx0)            (vx1)
   *     --o---------o------o------
   *       |         |      |      x
   *       0         p      1
   *
   */
  static T Linear(const T px, const T vx[2]);

  /**
   * Bilinear interpolation at point p
   *
   *     y
   *       |
   *       |
   *    1--o----------------o
   *       |(vy0)    |      |(vy1)
   *       |         |      |
   *   py--|---------p------|
   *       |         |      |
   *       |         |      |
   *       |         |      |
   *    0--o----------------o------
   *       |(vx0)    |      |(vx1) x
   *       0         px     1
   *
   */
  static T Bilinear(const T px, const T py, const T vx[2], const T vy[2]);

  /**
   * Trilinear interpolation at point p
   *
   *            o----------------o
   *           /|(vz2)          /|(vz3)
   *     y    / |              / |
   *       | /  |             /  |
   *       |/   |            /   |
   *    1--o----------------o    |
   *  (vy0)|    |/     (vy1)|    |
   *       | 1--o-----------|----o
   *       |   / (vz0)      |   / (vz1)
   *       |  /             |  /
   *       | /              | /
   *       |/               |/
   *    0--o----------------o------
   *      /|(vx0)           |(vx1) x
   *   z / 0                1
   *
   */
  static T Trilinear(const T px, const T py, const T pz, const T vx[2], const T vy[2], const T vz[4]);
};


} }
