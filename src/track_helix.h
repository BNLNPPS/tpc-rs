#ifndef TPCRS_TRACK_HELIX_H_
#define TPCRS_TRACK_HELIX_H_

#include <cmath>
#include <ostream>
#include <utility>
#include <algorithm>

#include "particles/StThreeVector.hh"


class TrackHelix
{
 public:

  TrackHelix() = default;

  /// curvature, dip angle, phase, origin, h
  TrackHelix(double c, double dip, double phase, const StThreeVector<double> &o, int h = -1);

  /// B: Signed magnetic field in kilogausses
  /// q: Charge of particle (+/- 1)
  TrackHelix(const StThreeVector<double> &p, const StThreeVector<double> &o, double B, double q);

  StThreeVector<double> momentum(double) const;     // returns the momentum at origin
  StThreeVector<double> momentumAt(double, double) const; // returns momemtum at S

  /// Returns charge of particle
  int charge(double B) const { return (B > 0 ? -h_ : h_); }

  // 2d DCA to x,y point signed relative to curvature
  double curvatureSignedDistance(double x, double y) ;
  // 2d DCA to x,y point signed relative to rotation
  double geometricSignedDistance(double x, double y) ;
  // 3d DCA to 3d point signed relative to curvature
  double curvatureSignedDistance(const StThreeVector<double> &) ;
  // 3d DCA to 3d point signed relative to rotation
  double geometricSignedDistance(const StThreeVector<double> &) ;

  double dipAngle() const { return dip_angle_; }

  /// 1/R in xy-plane
  double curvature() const { return curvature_; }

  /// aziumth in xy-plane measured from ring center
  double phase() const { return phase_; }

  /// x-center of circle in xy-plane
  double xcenter() const { return singularity_ ? 0 : origin_.x - cos_phase_ / curvature_; }

  /// y-center of circle in xy-plane
  double ycenter() const { return singularity_ ? 0 : origin_.y - sin_phase_ / curvature_; }

  /// -sign(q*B);
  int h() const { return h_; }

  /// Starting point
  const StThreeVector<double> &origin() const { return origin_; }

  void setParameters(double c, double dip, double phase, const StThreeVector<double> &o, int h);

  /// coordinates of helix at point s
  double x(double s) const
  {
    return singularity_ ?
      origin_.x - s * cos_dip_angle_ * sin_phase_ :
      origin_.x + (std::cos(phase_ + s * h_ * curvature_ * cos_dip_angle_) - cos_phase_) / curvature_;
  }

  double y(double s) const
  {
    return singularity_ ?
      origin_.y + s * cos_dip_angle_ * cos_phase_ :
      origin_.y + (std::sin(phase_ + s * h_ * curvature_ * cos_dip_angle_) - sin_phase_) / curvature_;
  }

  double z(double s) const { return origin_.z + s * sin_dip_angle_; }

  StThreeVector<double> at(double s) const { return StThreeVector<double>{x(s), y(s), z(s)}; }

  /// pointing vector of helix at point s
  double cx(double s) const
  {
    return singularity_ ?
      -cos_dip_angle_ * sin_phase_ :
      -std::sin(phase_ + s * h_ * curvature_ * cos_dip_angle_) * h_ * cos_dip_angle_;
  }

  double cy(double s) const
  {
    return singularity_ ?
      cos_dip_angle_ * cos_phase_ :
      std::cos(phase_ + s * h_ * curvature_ * cos_dip_angle_) * h_ * cos_dip_angle_;
  }

  double cz() const { return sin_dip_angle_; }

  StThreeVector<double> cat(double s) const { return StThreeVector<double>{cx(s), cy(s), cz()}; }

  /// returns period length of helix
  double period() const
  {
    return singularity_ ?
      std::numeric_limits<double>::max() :
      std::abs(2 * M_PI / (h_ * curvature_ * cos_dip_angle_));
  }

  /// path length at given r (cylindrical r)
  std::pair<double, double> pathLength(double r)   const;

  /// path length at given r (cylindrical r, cylinder axis at x,y)
  std::pair<double, double> pathLength(double r, double x, double y);

  /// path length at distance of closest approach to a given point
  double pathLength(const StThreeVector<double> &p, bool scanPeriods = true) const;

  /// path length at intersection with plane
  double pathLength(const StThreeVector<double> &r, const StThreeVector<double> &n) const;

  /// path length at distance of closest approach in the xy-plane to a given point
  double pathLength(double x, double y) const { return fudgePathLength(StThreeVector<double>{x, y, 0}); }

  /// path lengths in centimeters at dca between two helices
  std::pair<double, double> pathLengths(const TrackHelix &,
                                        double minStepSize = 1.e-3,
                                        double minRange = 10) const;

  /// minimal distance between point and helix
  double       distance(const StThreeVector<double> &p, bool scanPeriods = true) const;

  /// move the origin along the helix to s which becomes then s=0
  void moveOrigin(double s);

  static const double NoSolution;

 private:

  void setCurvature(double);
  void setPhase(double);
  void setDipAngle(double);

  /// value of S where distance in x-y plane is minimal
  double fudgePathLength(const StThreeVector<double> &) const;

  /// A flag indicating degenerate case of a straight line (B=0)
  bool singularity_;
  /// Track origin in cm
  StThreeVector<double> origin_;
  double dip_angle_;
  /// Track curvature measured in 1/cm
  double curvature_;
  double phase_;
  /// -sign(q*B);
  int h_;

  double cos_dip_angle_;
  double sin_dip_angle_;
  double cos_phase_;
  double sin_phase_;
};

int operator== (const TrackHelix &, const TrackHelix &);
int operator!= (const TrackHelix &, const TrackHelix &);
std::ostream &operator<<(std::ostream &, const TrackHelix &);

#endif
