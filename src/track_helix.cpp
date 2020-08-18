#include <limits>

#include "tpcrs/detail/track_helix.h"

const double TrackHelix::NoSolution = 3.e+33;

static const double c_light = 2.99792458e+10; // meter/second;

TrackHelix::TrackHelix(double c, double d, double phase, const Coords &o, int h)
{
  setParameters(c, d, phase, o, h);
}


/**
 * Momentum is in GeV
 * Origin is in cm
 * Magnetic field value B must be in kilogausses.
 * Charge of particle in positron charges (e.g. -1 for electron)
 */
TrackHelix::TrackHelix(const Coords &p, const Coords &o, double B, double q)
{
  h_ = (q * B <= 0) ? 1 : -1;

  if (p.y == 0 && p.x == 0)
    setPhase((M_PI_4) * (1 - 2.*h_));
  else
    setPhase(atan2(p.y, p.x) - h_ * M_PI_2);

  setDipAngle(atan2(p.z, p.perp()));
  origin_ = o;

  setCurvature(std::abs(c_light * q * B / p.mag() / cos_dip_angle_));
}


Coords TrackHelix::momentum(double B) const
{
  if (singularity_)
    return (Coords{0, 0, 0});
  else {
    {
      double pt = std::abs(c_light * B / curvature_);

      return (Coords{pt * std::cos(phase_ + h_* M_PI_2), // pos part pos field
                     pt * std::sin(phase_ + h_* M_PI_2),
                     pt * std::tan(dip_angle_)});
    }
  }
}


Coords TrackHelix::momentumAt(double S, double B) const
{
  // Obtain phase-shifted momentum from phase-shift of origin
  TrackHelix tmp(*this);
  tmp.moveOrigin(S);
  return tmp.momentum(B);
}


double TrackHelix::geometricSignedDistance(double x, double y)
{
  // Geometric signed distance
  double thePath = pathLength(x, y);
  Coords DCA2dPosition = at(thePath);
  DCA2dPosition.z = 0;
  Coords position{x, y, 0};
  Coords DCAVec = (DCA2dPosition - position);
  Coords momVec;

  // Deal with straight tracks
  if (singularity_) {
    momVec = at(1) - at(0);
    momVec.z = 0;
  }
  else {
    momVec = momentumAt(thePath, 1e13); // Don't care about Bmag.  Helicity is what matters.
    momVec.z = 0;
  }

  double cross = DCAVec.x * momVec.y - DCAVec.y * momVec.x;
  double theSign = (cross >= 0) ? 1. : -1.;
  return theSign * DCAVec.perp();
}


double TrackHelix::curvatureSignedDistance(double x, double y)
{
  // Protect against h_ = 0 or zero field
  if (singularity_ || std::abs(h_) <= 0) {
    return (geometricSignedDistance(x, y));
  }
  else {
    return (geometricSignedDistance(x, y)) / (h_);
  }
}


double TrackHelix::geometricSignedDistance(const Coords &pos)
{
  double sdca2d = geometricSignedDistance(pos.x, pos.y);
  double theSign = (sdca2d >= 0) ? 1. : -1.;
  return (distance(pos)) * theSign;
}


double TrackHelix::curvatureSignedDistance(const Coords &pos)
{
  double sdca2d = curvatureSignedDistance(pos.x, pos.y);
  double theSign = (sdca2d >= 0) ? 1. : -1.;
  return (distance(pos)) * theSign;
}


void TrackHelix::setParameters(double c, double dip, double phase,
                            const Coords &o, int h)
{
  // The order in which the parameters are set is important since setCurvature
  // might have to adjust the others.
  h_ = (h >= 0) ? 1 : -1;  // Default is: positive particle

  // positive field
  origin_   = o;
  setDipAngle(dip);
  setPhase(phase);

  // Check for singularity and correct for negative curvature. May change h_ and
  // phase_. Must therefore be set last.
  setCurvature(c);

  // For the case B=0, h is ill defined. In the following we always assume
  // h = +1. Since phase = psi - h * pi/2 we have to correct the phase in case
  // h = -1. This assumes that the user uses the same h for phase as the one he
  // passed to the constructor.
  if (singularity_ && h_ == -1) {
    h_ = +1;
    setPhase(phase_ - M_PI);
  }
}


void TrackHelix::setCurvature(double val)
{
  if (val < 0) {
    curvature_ = -val;
    h_ = -h_;
    setPhase(phase_ + M_PI);
  }
  else
    curvature_ = val;

  if (std::abs(curvature_) <= std::numeric_limits<double>::epsilon())
    singularity_ = true;   // straight line
  else
    singularity_ = false;  // curved
}


void TrackHelix::setPhase(double val)
{
  phase_       = val;
  cos_phase_    = std::cos(phase_);
  sin_phase_    = std::sin(phase_);

  if (std::abs(phase_) > M_PI)
    phase_ = std::atan2(sin_phase_, cos_phase_);  // force range [-pi,pi]
}


void TrackHelix::setDipAngle(double val)
{
  dip_angle_    = val;
  cos_dip_angle_ = std::cos(dip_angle_);
  sin_dip_angle_ = std::sin(dip_angle_);
}


double TrackHelix::fudgePathLength(const Coords &p) const
{
  double dx = p.x - origin_.x;
  double dy = p.y - origin_.y;

  double s;

  if (singularity_) {
    s = (dy * cos_phase_ - dx * sin_phase_) / cos_dip_angle_;
  }
  else {
    s = std::atan2(        dy * cos_phase_ - dx * sin_phase_,
          1 / curvature_ + dx * cos_phase_ + dy * sin_phase_) / (h_ * curvature_ * cos_dip_angle_);
  }

  return s;
}


double TrackHelix::distance(const Coords &p, bool scanPeriods) const
{
  return (at(pathLength(p, scanPeriods)) - p).mag();
}


/**
 * Calculates the path length at the distance of closest approach between the
 * helix and point p. For the case of B=0 (straight line) the path length can be
 * calculated analytically. For B>0 there is unfortunately no easy solution to
 * the problem. The Newton method is used to find the root of the referring
 * equation. The 'fudgePathLength' serves as a starting value.
 */
double TrackHelix::pathLength(const Coords &p, bool scanPeriods) const
{
  double s;
  double dx = p.x - origin_.x;
  double dy = p.y - origin_.y;
  double dz = p.z - origin_.z;

  if (singularity_) {
    s = cos_dip_angle_ * (cos_phase_ * dy - sin_phase_ * dx) +
        sin_dip_angle_ * dz;
  }
  else {
    const double MaxPrecisionNeeded = 1e-4; // micrometer
    const int    MaxIterations      = 100;

    // The math is taken from Maple with C(expr,optimized) and some
    // hand-editing. It is not very nice but efficient.
    double t34 = curvature_ * cos_dip_angle_ * cos_dip_angle_;
    double t41 = sin_dip_angle_ * sin_dip_angle_;
    double t6, t7, t11, t12, t19;

    // Get a first guess by using the dca in 2D. Since in some extreme cases we
    // might be off by n periods we add (subtract) periods in case we get any
    // closer.
    s = fudgePathLength(p);

    if (scanPeriods) {
      double period_len = period();
      int    j, jmin = 0;
      double d, dmin = (at(s) - p).mag();

      for (j = 1; j < MaxIterations; j++) {
        if ((d = (at(s + j * period_len) - p).mag()) < dmin) {
          dmin = d;
          jmin = j;
        }
        else
          break;
      }

      for (j = -1; -j < MaxIterations; j--) {
        if ((d = (at(s + j * period_len) - p).mag()) < dmin) {
          dmin = d;
          jmin = j;
        }
        else
          break;
      }

      if (jmin) s += jmin * period_len;
    }

    // Newtons method:
    // Stops after MaxIterations iterations or if the required precision is
    // obtained. Whatever comes first.
    double sOld = s;

    for (int i = 0; i < MaxIterations; i++) {
      t6  = phase_ + s * h_ * curvature_ * cos_dip_angle_;
      t7  = std::cos(t6);
      t11 = dx - (1 / curvature_) * (t7 - cos_phase_);
      t12 = std::sin(t6);
      t19 = dy - (1 / curvature_) * (t12 - sin_phase_);
      s  -= (t11 * t12 * h_ * cos_dip_angle_ - t19 * t7 * h_ * cos_dip_angle_ -
             (dz - s * sin_dip_angle_) * sin_dip_angle_) /
            (t12 * t12 * cos_dip_angle_ * cos_dip_angle_ + t11 * t7 * t34 +
             t7 * t7 * cos_dip_angle_ * cos_dip_angle_ +
             t19 * t12 * t34 + t41);

      if (std::abs(sOld - s) < MaxPrecisionNeeded) break;

      sOld = s;
    }
  }

  return s;
}


std::pair<double, double> TrackHelix::pathLength(double r) const
{
  std::pair<double, double> value;
  std::pair<double, double> VALUE(999999999., 999999999.);

  //
  // The math is taken from Maple with C(expr,optimized) and
  // some hand-editing. It is not very nice but efficient.
  // 'first' is the smallest of the two solutions (may be negative)
  // 'second' is the other.
  //
  if (singularity_) {
    double t1 = cos_dip_angle_ * (origin_.x * sin_phase_ - origin_.y * cos_phase_);
    double t12 = origin_.y * origin_.y;
    double t13 = cos_phase_ * cos_phase_;
    double t15 = r * r;
    double t16 = origin_.x * origin_.x;
    double t20 = -cos_dip_angle_ * cos_dip_angle_ * (2.0 * origin_.x * sin_phase_ * origin_.y * cos_phase_ +
                 t12 - t12 * t13 - t15 + t13 * t16);

    if (t20 < 0.) return VALUE;

    t20 = ::sqrt(t20);
    value.first  = (t1 - t20) / (cos_dip_angle_ * cos_dip_angle_);
    value.second = (t1 + t20) / (cos_dip_angle_ * cos_dip_angle_);
  }
  else {
    double t1 = origin_.y * curvature_;
    double t2 = sin_phase_;
    double t3 = curvature_ * curvature_;
    double t4 = origin_.y * t2;
    double t5 = cos_phase_;
    double t6 = origin_.x * t5;
    double t8 = origin_.x * origin_.x;
    double t11 = origin_.y * origin_.y;
    double t14 = r * r;
    double t15 = t14 * curvature_;
    double t17 = t8 * t8;
    double t19 = t11 * t11;
    double t21 = t11 * t3;
    double t23 = t5 * t5;
    double t32 = t14 * t14;
    double t35 = t14 * t3;
    double t38 = 8.0 * t4 * t6 - 4.0 * t1 * t2 * t8 - 4.0 * t11 * curvature_ * t6 +
                 4.0 * t15 * t6 + t17 * t3 + t19 * t3 + 2.0 * t21 * t8 + 4.0 * t8 * t23 -
                 4.0 * t8 * origin_.x * curvature_ * t5 - 4.0 * t11 * t23 -
                 4.0 * t11 * origin_.y * curvature_ * t2 + 4.0 * t11 - 4.0 * t14 +
                 t32 * t3 + 4.0 * t15 * t4 - 2.0 * t35 * t11 - 2.0 * t35 * t8;
    double t40 = (-t3 * t38);

    if (t40 < 0.) return VALUE;

    t40 = ::sqrt(t40);

    double t43 = origin_.x * curvature_;
    double t45 = 2.0 * t5 - t35 + t21 + 2.0 - 2.0 * t1 * t2 - 2.0 * t43 - 2.0 * t43 * t5 + t8 * t3;
    double t46 = h_ * cos_dip_angle_ * curvature_;

    value.first = (-phase_ + 2.0 * std::atan((-2.0 * t1 + 2.0 * t2 + t40) / t45)) / t46;
    value.second = -(phase_ + 2.0 * std::atan((2.0 * t1 - 2.0 * t2 + t40) / t45)) / t46;

    //
    //   Solution can be off by +/- one period, select smallest
    //
    double p = period();

    if (! std::isnan(value.first)) {
      if (std::abs(value.first - p) < std::abs(value.first)) value.first = value.first - p;
      else if (std::abs(value.first + p) < std::abs(value.first)) value.first = value.first + p;
    }

    if (! std::isnan(value.second)) {
      if (std::abs(value.second - p) < std::abs(value.second)) value.second = value.second - p;
      else if (std::abs(value.second + p) < std::abs(value.second)) value.second = value.second + p;
    }
  }

  if (value.first > value.second)
    std::swap(value.first, value.second);

  return (value);
}


std::pair<double, double> TrackHelix::pathLength(double r, double x, double y)
{
  double x0 = origin_.x;
  double y0 = origin_.y;
  origin_.x = x0 - x;
  origin_.y = y0 - y;
  std::pair<double, double> result = pathLength(r);
  origin_.x = x0;
  origin_.y = y0;
  return result;
}


/**
 * Vector 'r' defines the position of the center and vector 'n' the normal
 * vector of the plane.  For a straight line there is a simple analytical
 * solution. For curvatures > 0 the root is determined by Newton method. In case
 * no valid s can be found the max. largest value for s is returned.
 */
double TrackHelix::pathLength(const Coords &r, const Coords &n) const
{
  double s;

  if (singularity_) {
    double t = n.z * sin_dip_angle_ +
               n.y * cos_dip_angle_ * cos_phase_ -
               n.x * cos_dip_angle_ * sin_phase_;

    if (t == 0)
      s = NoSolution;
    else
      s = ((r - origin_) * n) / t;
  }
  else {
    const double MaxPrecisionNeeded = 1e-4; // micrometer;
    const int    MaxIterations      = 20;

    double A = curvature_ * ((origin_ - r) * n) -
               n.x * cos_phase_ -
               n.y * sin_phase_;
    double t = h_ * curvature_ * cos_dip_angle_;
    double u = n.z * curvature_ * sin_dip_angle_;

    double a, f, fp;
    double sOld = s = 0;
    double shiftOld = 0;
    double shift;
    //		(std::cos(angMax)-1)/angMax = 0.1
    const double angMax = 0.21;
    double deltas = std::abs(angMax / (curvature_ * cos_dip_angle_));
    //              dampingFactor = exp(-0.5);
    //	double dampingFactor = 0.60653;
    int i;

    for (i = 0; i < MaxIterations; i++) {
      a  = t * s + phase_;
      double sina = std::sin(a);
      double cosa = std::cos(a);
      f  = A +
           n.x * cosa +
           n.y * sina +
           u * s;
      fp = -n.x * sina * t +
           n.y * cosa * t +
           u;

      if ( std::abs(fp)*deltas <= std::abs(f) ) { //too big step
        int sgn = 1;

        if (fp < 0.) sgn = -sgn;

        if (f < 0.) sgn = -sgn;

        shift = sgn * deltas;

        if (shift < 0) shift *= 0.9; // don't get stuck shifting +/-deltas
      }
      else {
        shift = f / fp;
      }

      s -= shift;
      shiftOld = shift;

      if (std::abs(sOld - s) < MaxPrecisionNeeded) break;

      sOld = s;
    }

    if (i == MaxIterations) return NoSolution;
  }

  return s;
}


std::pair<double, double>
TrackHelix::pathLengths(const TrackHelix &h, double minStepSize, double minRange) const
{
  // Cannot handle case where one is a helix and the other one is a straight line
  if (singularity_ != h.singularity_)
    return std::pair<double, double>(NoSolution, NoSolution);

  double s1, s2;

  if (singularity_) {
    //
    //  Analytic solution
    //
    Coords dv = h.origin_ - origin_;
    Coords a{-cos_dip_angle_ * sin_phase_,
                            cos_dip_angle_ * cos_phase_,
                            sin_dip_angle_};
    Coords b{-h.cos_dip_angle_ * h.sin_phase_,
                            h.cos_dip_angle_ * h.cos_phase_,
                            h.sin_dip_angle_};
    double ab = a * b;
    double g  = dv * a;
    double k  = dv * b;
    s2 = (k - ab * g) / (ab * ab - 1.);
    s1 = g + s2 * ab;
    return std::pair<double, double>(s1, s2);
  }
  else {
    // First step: get dca in the xy-plane as start value
    double dx = h.xcenter() - xcenter();
    double dy = h.ycenter() - ycenter();
    double dd = ::sqrt(dx * dx + dy * dy);
    double r1 = 1 / curvature();
    double r2 = 1 / h.curvature();

    double cosAlpha = (r1 * r1 + dd * dd - r2 * r2) / (2 * r1 * dd);

    double s;
    double x, y;

    if (std::abs(cosAlpha) < 1) {           // two solutions
      double sinAlpha = std::sin(std::acos(cosAlpha));
      x = xcenter() + r1 * (cosAlpha * dx - sinAlpha * dy) / dd;
      y = ycenter() + r1 * (sinAlpha * dx + cosAlpha * dy) / dd;
      s = pathLength(x, y);
      x = xcenter() + r1 * (cosAlpha * dx + sinAlpha * dy) / dd;
      y = ycenter() + r1 * (cosAlpha * dy - sinAlpha * dx) / dd;
      double a = pathLength(x, y);

      if (h.distance(at(a)) < h.distance(at(s))) s = a;
    }
    else {                              // no intersection (or exactly one)
      int rsign = ((r2 - r1) > dd ? -1 : 1); // set -1 when *this* helix is
      // completely contained in the other
      x = xcenter() + rsign * r1 * dx / dd;
      y = ycenter() + rsign * r1 * dy / dd;
      s = pathLength(x, y);
    }

    // Second step: scan in decreasing intervals around seed 's' minRange and
    // minStepSize are passed as arguments to the method. They have default
    // values defined in the header file.
    double dmin  = h.distance(at(s));
    double range = std::max(2 * dmin, minRange);
    double ds    = range / 10;
    double slast = -999999, ss, d;
    s1 = s - range / 2.;
    s2 = s + range / 2.;

    while (ds > minStepSize) {
      for (ss = s1; ss < s2 + ds; ss += ds) {
        d = h.distance(at(ss));

        if (d < dmin) {
          dmin = d;
          s = ss;
        }

        slast = ss;
      }

      // In the rare cases where the minimum is at the the border of the current
      // range we shift the range and start all over, i.e we do not decrease
      // 'ds'.  Else we decrease the search intervall around the current minimum
      // and redo the scan in smaller steps.
      if (s == s1) {
        d = 0.8 * (s2 - s1);
        s1 -= d;
        s2 -= d;
      }
      else if (s == slast) {
        d = 0.8 * (s2 - s1);
        s1 += d;
        s2 += d;
      }
      else {
        s1 = s - ds;
        s2 = s + ds;
        ds /= 10;
      }
    }

    return std::pair<double, double>(s, h.pathLength(at(s)));
  }
}


void TrackHelix::moveOrigin(double s)
{
  if (singularity_)
    origin_	= at(s);
  else {
    Coords newOrigin = at(s);
    double newPhase = std::atan2(newOrigin.y - ycenter(),
                            newOrigin.x - xcenter());
    origin_ = newOrigin;
    setPhase(newPhase);
  }
}


int operator== (const TrackHelix &a, const TrackHelix &b)
{
  return (a.origin()    == b.origin()    &&
          a.dipAngle()  == b.dipAngle()  &&
          a.curvature() == b.curvature() &&
          a.phase()     == b.phase()     &&
          a.h()         == b.h());
}

int operator!= (const TrackHelix &a, const TrackHelix &b) {return !(a == b);}

std::ostream &operator<<(std::ostream &os, const TrackHelix &h)
{
  return os << '('
         << "curvature = "  << h.curvature() << ", "
         << "dip angle = "  << h.dipAngle()  << ", "
         << "phase = "      << h.phase()     << ", "
         << "h = "          << h.h()         << ", "
         << "origin = "     << h.origin().x << " " << h.origin().y << " " << h.origin().z << ')';
}
