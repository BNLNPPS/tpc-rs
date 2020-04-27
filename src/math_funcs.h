#ifndef TPCRS_MATH_H_
#define TPCRS_MATH_H_

namespace tpcrs {

/// Round to nearest integer. Rounds half integers to the nearest even integer.
///
/// This code is copied from the ROOT project v6.20.04 covered by the LGPL
/// https://root.cern/  https://github.com/root-project/root
/// See TMath::Nint(T x) in root/math/mathcore/inc/TMath.h
template<typename T>
int irint(T x)
{
   int i;
   if (x >= 0) {
      i = int(x + 0.5);
      if ( i & 1 && x + 0.5 == T(i) ) i--;
   } else {
      i = int(x - 0.5);
      if ( i & 1 && x - 0.5 == T(i) ) i++;
   }
   return i;
}

// Bessel functions
double BesselI0(double x);         /// modified Bessel function I_0(x)
double BesselK0(double x);         /// modified Bessel function K_0(x)
double BesselI1(double x);         /// modified Bessel function I_1(x)
double BesselK1(double x);         /// modified Bessel function K_1(x)

double Gamma(double a,double x);
double GammaDist(double x, double gamma, double mu=0, double beta=1);

}

#endif
