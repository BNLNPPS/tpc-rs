/***************************************************************************
 * Author: Brian Lasiuk, Sep 1997
 *
 * Parametrization of a physical helix. See the SCL user guide for more.
 **************************************************************************/
#ifndef ST_PHYSICAL_HELIX_HH
#define ST_PHYSICAL_HELIX_HH

#include "StarClassLibrary/StThreeVector.hh"
#include "StarClassLibrary/StHelix.hh"

class StPhysicalHelix : public StHelix
{
 public:
  // Requires: momentum, origin, signed Magnetic Field in kilogausses
  //           and Charge of particle (+/- 1)
  StPhysicalHelix(const StThreeVector<double> &,
                  const StThreeVector<double> &,
                  double, double);

  // curvature, dip angle, phase, origin, h
  StPhysicalHelix(double, double, double,
                  const StThreeVector<double> &, int h = -1);
  StPhysicalHelix();

  ~StPhysicalHelix();

  // Requires:  signed Magnetic Field
  StThreeVector<double> momentum(double) const;     // returns the momentum at origin
  StThreeVector<double> momentumAt(double, double) const; // returns momemtum at S
  int                   charge(double)   const;     // returns charge of particle
  // 2d DCA to x,y point signed relative to curvature
  double curvatureSignedDistance(double x, double y) ;
  // 2d DCA to x,y point signed relative to rotation
  double geometricSignedDistance(double x, double y) ;
  // 3d DCA to 3d point signed relative to curvature
  double curvatureSignedDistance(const StThreeVector<double> &) ;
  // 3d DCA to 3d point signed relative to rotation
  double geometricSignedDistance(const StThreeVector<double> &) ;
};

#endif
