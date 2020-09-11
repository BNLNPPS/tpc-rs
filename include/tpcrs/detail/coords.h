#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "TGeoMatrix.h"

#include "tpcrs/detail/struct_containers.h"


template<class T> struct ThreeVector
{
  T x, y, z;

  void setPhi(T angle)
  {
    double  r = mag();
    double th = theta();

    x = r * std::sin(th) * std::cos(angle);
    y = r * std::sin(th) * std::sin(angle);
  }

  void setTheta(T angle)
  {
    double r  = mag();
    double ph = phi();

    x = r * std::sin(angle) * std::cos(ph);
    y = r * std::sin(angle) * std::sin(ph);
    z = r * std::cos(angle);
  }

  void setMagnitude(T r)
  {
    double th = theta();
    double ph = phi();

    x = r * std::sin(th) * std::cos(ph);
    y = r * std::sin(th) * std::sin(ph);
    z = r * std::cos(th);
  }

  const T* xyz() const { return &x; }
  T* xyz() { return &x; }
  T  theta() const { return std::acos(cosTheta()); }
  T  cosTheta() const { return z / (mag() + 1e-20); }
  T  phi() const { return std::atan2(y, x); }
  T  perp() const { return std::sqrt(x * x + y * y); }
  T  perp2() const { return x * x + y * y; }
  T  mag() const { return std::sqrt(x * x + y * y + z * z); }
  T  mag2() const { return x * x + y * y + z * z; }

  T  pseudoRapidity() const
  {
    // return 0.5*::log( (mag() + z())/(mag() - z()) );
    double tmp = tan(theta() / 2.);
    if (tmp <= 0.) return 1e20;
    return -std::log(tmp);
  }

  T operator[] (size_t i) const
  {
    if (i <= 2)  return (&x)[i];
    throw std::out_of_range("ThreeVector<T>::operator[]: bad index");
    return T{};
  }

  T &operator[] (size_t i)
  {
    if (i <= 2)  return (&x)[i];
    throw std::out_of_range("ThreeVector<T>::operator[]: bad index");
    return x;
  }

  ThreeVector<T> unit() const
  {
    double tmp = mag();
    if (tmp <= 0.) tmp = 1e-20;
    return *this / tmp;
  }

  void rotateX(T angle)
  {
    double yPrime = std::cos(angle) * y - std::sin(angle) * z;
    double zPrime = std::sin(angle) * y + std::cos(angle) * z;

    y = yPrime;
    z = zPrime;
  }

  void rotateY(T angle)
  {
    double zPrime = std::cos(angle) * z - std::sin(angle) * x;
    double xPrime = std::sin(angle) * z + std::cos(angle) * x;

    x = xPrime;
    z = zPrime;
  }

  void rotateZ(T angle)
  {
    double xPrime = std::cos(angle) * x - std::sin(angle) * y;
    double yPrime = std::sin(angle) * x + std::cos(angle) * y;

    x = xPrime;
    y = yPrime;
  }

  ThreeVector<T> operator-() { return ThreeVector<T>(-x, -y, -z); }
  ThreeVector<T> operator+() { return *this; }

  ThreeVector<T> &operator*= (double c) { x *= c; y *= c; z *= c; return *this; }
  ThreeVector<T> &operator/= (double c) { x /= c; y /= c; z /= c; return *this; }

  T angle(const T &v) const
  {
    double norm = this->mag2() * v.mag2();
    return norm > 0 ? std::acos(this->dot(v) / (std::sqrt(norm))) : 0;
  }

  ThreeVector<T> cross(const T& v) const
  {
    return ThreeVector<T>(y * v.z - z * v.y,
                            z * v.x - x * v.z,
                            x * v.y - y * v.x);
  }

  T dot(const ThreeVector<T>& v) const { return x * v.x + y * v.y + z * v.z; }

  bool operator== (const ThreeVector<T> &v) const { return x == v.x && y == v.y && z == v.z; }
  bool operator!= (const ThreeVector<T> &v) const { return !(*this == v); }

  ThreeVector<T> &operator+=(const ThreeVector<T>& v)
  {
    x += v.x; y += v.y; z += v.z;
    return *this;
  }

  ThreeVector<T> &operator-=(const ThreeVector<T>& v)
  {
    x -= v.x; y -= v.y; z -= v.z;
    return *this;
  }
};


// Non-member functions
template<class T>
inline ThreeVector<T>
operator+ (const ThreeVector<T> &v1, const ThreeVector<T> &v2)
{
  return ThreeVector<T>(v1) += v2;
}

template<class T>
inline ThreeVector<T>
operator-(const ThreeVector<T> &v1, const ThreeVector<T> &v2)
{
  return ThreeVector<T>(v1) -= v2;
}

template<class T>
inline T operator* (const ThreeVector<T> &v1, const ThreeVector<T> &v2)
{
  return ThreeVector<T>(v1).dot(v2);
}

template<class T>
inline ThreeVector<T> operator*(const ThreeVector<T> &v, double c)
{
  return ThreeVector<T>(v) *= c;
}

template<class T>
inline ThreeVector<T> operator*(double c, const ThreeVector<T> &v)
{
  return ThreeVector<T>(v) *= c;
}

template<class T>
inline ThreeVector<T> operator/(const ThreeVector<T> &v, double c)
{
  return ThreeVector<T>(v) /= c;
}

using Coords = ThreeVector<double>;

double mag(const Coords& c);
double perp(const Coords& c);
Coords unit(const Coords& c);
Coords operator-(const Coords &c1, const Coords &c2);
double operator*(const Coords &c1, const Coords &c2);
Coords operator/(const Coords &c, double v);

struct StTpcPadCoordinate
{
  int   sector;
  int   row;
  float pad;
  float timeBucket;
};

template<int dummy>
struct GCoordinate
{
  Coords position;
};

template<int dummy>
struct TCoordinate
{
  Coords position;
  int sector;
  int row;
};

using StGlobalCoordinate = GCoordinate<1>;
using StGlobalDirection = GCoordinate<2>;

using StTpcCoordinate = TCoordinate<1>;
using StTpcLocalCoordinate = TCoordinate<2>;
using StTpcLocalDirection = TCoordinate<3>;
using StTpcLocalSectorCoordinate = TCoordinate<4>;
using StTpcLocalSectorDirection = TCoordinate<5>;

// pad              => sector12       =>   subsector => sector => tpc      => global
// TpcPadCoordinate => TpcSectL => TpcSectLAligned => TpcLocal => Global
struct CoordTransform
{
  CoordTransform(const tpcrs::Configurator& cfg);

  // Raw Data <--> Tpc Local Sector Coordinates
  void local_sector_to_hardware(const StTpcLocalSectorCoordinate &a, StTpcPadCoordinate &b, bool useT0 = false, bool useTau = true) const;
  void hardware_to_local_sector(const StTpcPadCoordinate &a, StTpcLocalSectorCoordinate &b, bool useT0 = false, bool useTau = true) const;

  // Raw Data <--> Tpc Local Coordinates
  void local_to_hardware(const StTpcLocalCoordinate &a, StTpcPadCoordinate &b, bool useT0 = false, bool useTau = true) const
  {
    StTpcLocalSectorCoordinate c;
    local_to_local_sector(a, c);
    local_sector_to_hardware(c, b, useT0, useTau);
  }

  void hardware_to_local(const StTpcPadCoordinate &a, StTpcLocalCoordinate &b, bool useT0 = false, bool useTau = true) const
  {
    StTpcLocalSectorCoordinate c;
    hardware_to_local_sector(a, c, useT0, useTau);
    local_sector_to_local(c, b);
  }

  // Tpc Local Sector <--> TPC Local
  void local_sector_to_local(const StTpcLocalSectorCoordinate &a, StTpcLocalCoordinate &b) const;

  void local_sector_to_local_dir(const StTpcLocalSectorDirection &a, StTpcLocalDirection &b) const
  {
    Pad2Tpc(a.sector, a.row).LocalToMasterVect(a.position.xyz(), b.position.xyz());
    b.sector = a.sector;
    b.row = a.row;
  }

  void local_sector_to_global(const StTpcLocalSectorCoordinate &a, StGlobalCoordinate &b) const
  {
    StTpcLocalCoordinate c;
    local_sector_to_local(a, c);
    local_to_global(c, b);
  }

  void local_sector_to_global_dir(const StTpcLocalSectorDirection &a, StGlobalDirection &b) const
  {
    Pad2Glob(a.sector, a.row).LocalToMasterVect(a.position.xyz(), b.position.xyz());
  }

  void local_to_local_sector(const StTpcLocalCoordinate &a, StTpcLocalSectorCoordinate &b) const;

  void local_to_local_sector_dir(const StTpcLocalDirection &a, StTpcLocalSectorDirection &b) const
  {
    Pad2Tpc(a.sector, a.row).MasterToLocalVect(a.position.xyz(), b.position.xyz());
    b.sector = a.sector;
    b.row = a.row;
  }

  void global_to_local_sector(const StGlobalCoordinate &a, StTpcLocalSectorCoordinate &b, int sector, int row) const
  {
    StTpcLocalCoordinate c;
    global_to_local(a, c, sector, row);
    local_to_local_sector(c, b);
  }

  void global_to_local_sector_dir(const  StGlobalDirection &a, StTpcLocalSectorDirection &b, int sector, int row) const
  {
    Pad2Glob(sector, row).MasterToLocalVect(a.position.xyz(), b.position.xyz());
    b.sector = sector;
    b.row = row;
  }

  // Internal TpcCoordinate <-->  Global Coordinate
  void local_to_global(const StTpcLocalCoordinate &a, StGlobalCoordinate &b) const
  {
    tpc2global_.LocalToMaster(a.position.xyz(), b.position.xyz());
  }

  void global_to_local(const StGlobalCoordinate &a, StTpcLocalCoordinate &b, int sector, int row) const
  {
    tpc2global_.MasterToLocal(a.position.xyz(), b.position.xyz());
    b.sector = sector;
    b.row = row;
  }

  void local_to_global_dir(const StTpcLocalDirection &a, StGlobalDirection &b) const
  {
    tpc2global_.LocalToMasterVect(a.position.xyz(), b.position.xyz());
  }

  void global_to_local_dir(const StGlobalDirection &a, StTpcLocalDirection &b, int sector, int row) const
  {
    tpc2global_.MasterToLocalVect(a.position.xyz(), b.position.xyz());
    b.sector = sector;
    b.row = row;
  }

  // Raw Data <-->  Global Coordinate
  void hardware_to_global(const StTpcPadCoordinate &a, StGlobalCoordinate &b, bool useT0 = false, bool useTau = true) const
  {
    StTpcLocalCoordinate c;
    hardware_to_local(a, c, useT0, useTau);
    local_to_global(c, b);
  }

  void global_to_hardware(const StGlobalCoordinate &a, StTpcPadCoordinate &b, int sector, int row, bool useT0 = false, bool useTau = true) const
  {
    StTpcLocalCoordinate c;
    global_to_local(a, c, sector, row);
    local_to_hardware(c, b, useT0, useTau);
  }

  double ZToTime(double z, int sector, int row, int pad = 0) const;
  double TimeToZ(double tb, int sector, int row, int pad = 0) const;

  // Raw Data (pad row timebin or drift L From tpc local sector Coordinates
  int YToRow(double y, int sector) const;

  double XToPad(double x, int sector, int row) const;
  double PadToX(int sector, int row, double pad) const;

 private:

  // Glob     = Global coordinate
  // Tpc      = Tpc    -"-                survey
  // Half     = Tpc Half west / east -"-  survey
  // SupS     = super sector misalignment(?)
  // SubS[io] = SubSector[io] misalignment
  // SecL     = sector -"- coordinate (y_p, x_p, DriftDistance - z_p);
  // Pad      = Pad -"- (x_p,y_p,z_p) (Sector12 coordinate system)
  // Tpc => Global is tpc2global_
  // Pad => SecL   is internal Flip matrix
  enum ETpcSectorRotationType {kUndefSector     = -2,
                               kFlip            = -1, // Flip * Subs[io] => SupS
                               kSupS2Tpc        = 0, // SupS => Tpc
                               kSupS2Glob       = 1, // SupS => Tpc => Glob;
                               kSubSInner2SupS  = 2, // Subs[io] => SupS
                               kSubSOuter2SupS  = 3, // -"-
                               kSubSInner2Tpc   = 4, // (Subs[io] => SupS) => Tpc
                               kSubSOuter2Tpc   = 5, // -"-
                               kSubSInner2Glob  = 6, // (Subs[io] => SupS => Tpc) => Glob
                               kSubSOuter2Glob  = 7, // -"-
                               kPadInner2SupS   = 8, // (Pad => SecL) => (SubS[io] => SupS)
                               kPadOuter2SupS   = 9, // -"-
                               kPadInner2Tpc    = 10, // (Pad => SecL) => (SubS[io] => SupS => Tpc)
                               kPadOuter2Tpc    = 11, // -"-
                               kPadInner2Glob   = 12, // (Pad => SecL) => (SubS[io] => SupS => Tpc => Glob)
                               kPadOuter2Glob   = 13, // -"-
                               kTotalTpcSectorRotaions = 14
                              };

  const tpcrs::Configurator& cfg_;
  double timebin_width_;
  double z_inner_offset_;
  double z_outer_offset_;

  TGeoHMatrix tpc2global_;
  std::vector<TGeoHMatrix> sector_rotations_;

  void SetTpcRotations();

  const TGeoHMatrix &TpcRot(int sector, int k)      const {return sector_rotations_[kTotalTpcSectorRotaions*(sector - 1) + k];}
  const TGeoHMatrix &SupS2Tpc(int sector = 1)       const {return TpcRot(sector, kSupS2Tpc);}
  const TGeoHMatrix &SupS2Glob(int sector = 1)      const {return TpcRot(sector, kSupS2Glob);}
  const TGeoHMatrix &SubSInner2SupS(int sector = 1) const {return TpcRot(sector, kSubSInner2SupS);}

  const TGeoHMatrix &SubSOuter2SupS(int sector = 1) const {return TpcRot(sector, kSubSOuter2SupS);}
  const TGeoHMatrix &SubSInner2Tpc(int sector = 1)  const {return TpcRot(sector, kSubSInner2Tpc);}
  const TGeoHMatrix &SubSOuter2Tpc(int sector = 1)  const {return TpcRot(sector, kSubSOuter2Tpc);}
  const TGeoHMatrix &SubSInner2Glob(int sector = 1) const {return TpcRot(sector, kSubSInner2Glob);}
  const TGeoHMatrix &SubSOuter2Glob(int sector = 1) const {return TpcRot(sector, kSubSOuter2Glob);}

  const TGeoHMatrix &PadInner2SupS(int sector = 1)  const {return TpcRot(sector, kPadInner2SupS);}
  const TGeoHMatrix &PadOuter2SupS(int sector = 1)  const {return TpcRot(sector, kPadOuter2SupS);}
  const TGeoHMatrix &PadInner2Tpc(int sector = 1)   const {return TpcRot(sector, kPadInner2Tpc);}
  const TGeoHMatrix &PadOuter2Tpc(int sector = 1)   const {return TpcRot(sector, kPadOuter2Tpc);}
  const TGeoHMatrix &PadInner2Glob(int sector = 1)  const {return TpcRot(sector, kPadInner2Glob);}
  const TGeoHMatrix &PadOuter2Glob(int sector = 1)  const {return TpcRot(sector, kPadOuter2Glob);}

  const TGeoHMatrix &SubS2SupS(int sector = 1, int row = 1) const {return TpcRot(sector, tpcrs::IsInner(row, cfg_) ? kSubSInner2SupS : kSubSOuter2SupS);}
  const TGeoHMatrix &SubS2Tpc (int sector = 1, int row = 1) const {return TpcRot(sector, tpcrs::IsInner(row, cfg_) ? kSubSInner2Tpc  : kSubSOuter2Tpc );}
  const TGeoHMatrix &SubS2Glob(int sector = 1, int row = 1) const {return TpcRot(sector, tpcrs::IsInner(row, cfg_) ? kSubSInner2Glob : kSubSOuter2Glob);}

  const TGeoHMatrix &Pad2SupS(int sector = 1, int row = 1)  const {return TpcRot(sector, tpcrs::IsInner(row, cfg_) ? kPadInner2SupS : kPadOuter2SupS);}
  const TGeoHMatrix &Pad2Tpc (int sector = 1, int row = 1)  const {return TpcRot(sector, tpcrs::IsInner(row, cfg_) ? kPadInner2Tpc  : kPadOuter2Tpc );}
  const TGeoHMatrix &Pad2Glob(int sector = 1, int row = 1)  const {return TpcRot(sector, tpcrs::IsInner(row, cfg_) ? kPadInner2Glob : kPadOuter2Glob);}
};
