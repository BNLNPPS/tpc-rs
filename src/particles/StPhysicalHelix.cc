/***************************************************************************
 * Author: Brian Lasiuk, Sep 1997
 *
 * Parametrization of a physical helix.
 **************************************************************************/
#include <cmath>
#include "particles/StHelix.hh"
#include "particles/StPhysicalHelix.hh"
#include "particles/SystemOfUnits.h"

static const double c_light   = 2.99792458e+8 * meter/second;

StPhysicalHelix::StPhysicalHelix() {}

StPhysicalHelix::~StPhysicalHelix() { /* nop */ }

StPhysicalHelix::StPhysicalHelix(const StThreeVector<double> &p,
                                 const StThreeVector<double> &o,
                                 double B, double q)
{
  mH = (q * B <= 0) ? 1 : -1;

  if (p.y() == 0 && p.x() == 0)
    setPhase((M_PI / 4) * (1 - 2.*mH));
  else
    setPhase(atan2(p.y(), p.x()) - mH * M_PI / 2);

  setDipAngle(atan2(p.z(), p.perp()));
  mOrigin = o;

  {
    using namespace units;
    setCurvature(std::abs((c_light * nanosecond / meter * q * B / tesla) /
                      (abs(p) / GeV * mCosDipAngle) / meter));
  }
}

StPhysicalHelix::StPhysicalHelix(double c, double d, double phase,
                                 const StThreeVector<double> &o, int h)
  : StHelix(c, d, phase, o, h) { /* nop */}


StThreeVector<double> StPhysicalHelix::momentum(double B) const
{
  if (mSingularity)
    return (StThreeVector<double>(0, 0, 0));
  else {
    {
      using namespace units;
      double pt = GeV * std::abs(c_light* nanosecond / meter* B / tesla) / (std::abs(mCurvature) * meter);

      return (StThreeVector<double>(pt * cos(mPhase + mH* M_PI / 2), // pos part pos field
                                    pt * sin(mPhase + mH* M_PI / 2),
                                    pt * tan(mDipAngle)));
    }
  }
}

StThreeVector<double> StPhysicalHelix::momentumAt(double S, double B) const
{
  // Obtain phase-shifted momentum from phase-shift of origin
  StPhysicalHelix tmp(*this);
  tmp.moveOrigin(S);
  return tmp.momentum(B);
}

int StPhysicalHelix::charge(double B) const
{
  return (B > 0 ? -mH : mH);
}

double StPhysicalHelix::geometricSignedDistance(double x, double y)
{
  // Geometric signed distance
  double thePath = this->pathLength(x, y);
  StThreeVector<double> DCA2dPosition = this->at(thePath);
  DCA2dPosition.setZ(0);
  StThreeVector<double> position(x, y, 0);
  StThreeVector<double> DCAVec = (DCA2dPosition - position);
  StThreeVector<double> momVec;

  // Deal with straight tracks
  if (this->mSingularity) {
    momVec = this->at(1) - this->at(0);
    momVec.setZ(0);
  }
  else {
    momVec = this->momentumAt(thePath, 1. / tesla); // Don't care about Bmag.  Helicity is what matters.
    momVec.setZ(0);
  }

  double cross = DCAVec.x() * momVec.y() - DCAVec.y() * momVec.x();
  double theSign = (cross >= 0) ? 1. : -1.;
  return theSign * DCAVec.perp();
}

double StPhysicalHelix::curvatureSignedDistance(double x, double y)
{
  // Protect against mH = 0 or zero field
  if (this->mSingularity || abs(this->mH) <= 0) {
    return (this->geometricSignedDistance(x, y));
  }
  else {
    return (this->geometricSignedDistance(x, y)) / (this->mH);
  }

}

double StPhysicalHelix::geometricSignedDistance(const StThreeVector<double> &pos)
{
  double sdca2d = this->geometricSignedDistance(pos.x(), pos.y());
  double theSign = (sdca2d >= 0) ? 1. : -1.;
  return (this->distance(pos)) * theSign;
}

double StPhysicalHelix::curvatureSignedDistance(const StThreeVector<double> &pos)
{
  double sdca2d = this->curvatureSignedDistance(pos.x(), pos.y());
  double theSign = (sdca2d >= 0) ? 1. : -1.;
  return (this->distance(pos)) * theSign;
}
