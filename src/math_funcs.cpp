/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see https://root.cern/license                 *
 * For the list of contributors see root/README/CREDITS.                 *
 *************************************************************************/

// This code is copied from the ROOT project v6.20.04 covered by the LGPL
// 
// https://root.cern/
// https://github.com/root-project/root
//
// The original files:
//
// root/math/mathcore/src/TMath.cxx
// root/math/mathcore/inc/PdfFuncMathCore.h

#include <cmath>

#include "tpcrs/logger.h"
#include "math_cephes.h"

namespace tpcrs {

////////////////////////////////////////////////////////////////////////////////
/// Compute the modified Bessel function I_0(x) for any real x.
///
/// \author NvE 12-mar-2000 UU-SAP Utrecht
double BesselI0(double x)
{
   // Parameters of the polynomial approximation
   const double p1=1.0,          p2=3.5156229,    p3=3.0899424,
                  p4=1.2067492,    p5=0.2659732,    p6=3.60768e-2,  p7=4.5813e-3;

   const double q1= 0.39894228,  q2= 1.328592e-2, q3= 2.25319e-3,
                  q4=-1.57565e-3,  q5= 9.16281e-3,  q6=-2.057706e-2,
                  q7= 2.635537e-2, q8=-1.647633e-2, q9= 3.92377e-3;

   const double k1 = 3.75;
   double ax = std::abs(x);

   double y=0, result=0;

   if (ax < k1) {
      double xx = x/k1;
      y = xx*xx;
      result = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
   } else {
      y = k1/ax;
      result = (std::exp(ax)/std::sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the modified Bessel function K_0(x) for positive real x.
///
///  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
///     Applied Mathematics Series vol. 55 (1964), Washington.
///
/// \author NvE 12-mar-2000 UU-SAP Utrecht
double BesselK0(double x)
{
   // Parameters of the polynomial approximation
   const double p1=-0.57721566,  p2=0.42278420,   p3=0.23069756,
                  p4= 3.488590e-2, p5=2.62698e-3,   p6=1.0750e-4,    p7=7.4e-6;

   const double q1= 1.25331414,  q2=-7.832358e-2, q3= 2.189568e-2,
                  q4=-1.062446e-2, q5= 5.87872e-3,  q6=-2.51540e-3,  q7=5.3208e-4;

   if (x <= 0) {
      LOG_ERROR << "tpcrs::BesselK0: *K0* Invalid argument x = " << x << '\n';
      return 0;
   }

   double y=0, result=0;

   if (x <= 2) {
      y = x*x/4;
      result = (-std::log(x/2.)*BesselI0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = 2/x;
      result = (std::exp(-x)/std::sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the modified Bessel function I_1(x) for any real x.
///
///  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
///     Applied Mathematics Series vol. 55 (1964), Washington.
///
/// \author NvE 12-mar-2000 UU-SAP Utrecht
double BesselI1(double x)
{
   // Parameters of the polynomial approximation
   const double p1=0.5,          p2=0.87890594,   p3=0.51498869,
                  p4=0.15084934,   p5=2.658733e-2,  p6=3.01532e-3,  p7=3.2411e-4;

   const double q1= 0.39894228,  q2=-3.988024e-2, q3=-3.62018e-3,
                  q4= 1.63801e-3,  q5=-1.031555e-2, q6= 2.282967e-2,
                  q7=-2.895312e-2, q8= 1.787654e-2, q9=-4.20059e-3;

   const double k1 = 3.75;
   double ax = std::abs(x);

   double y=0, result=0;

   if (ax < k1) {
      double xx = x/k1;
      y = xx*xx;
      result = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = k1/ax;
      result = (std::exp(ax)/std::sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
      if (x < 0) result = -result;
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the modified Bessel function K_1(x) for positive real x.
///
///  M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
///     Applied Mathematics Series vol. 55 (1964), Washington.
///
/// \author NvE 12-mar-2000 UU-SAP Utrecht
double BesselK1(double x)
{
   // Parameters of the polynomial approximation
   const double p1= 1.,          p2= 0.15443144,  p3=-0.67278579,
                  p4=-0.18156897,  p5=-1.919402e-2, p6=-1.10404e-3,  p7=-4.686e-5;

   const double q1= 1.25331414,  q2= 0.23498619,  q3=-3.655620e-2,
                  q4= 1.504268e-2, q5=-7.80353e-3,  q6= 3.25614e-3,  q7=-6.8245e-4;

   if (x <= 0) {
      LOG_ERROR << "tpcrs::BesselK1: *K1* Invalid argument x = " << x << '\n';
      return 0;
   }

   double y=0,result=0;

   if (x <= 2) {
      y = x*x/4;
      result = (std::log(x/2.)*BesselI1(x))+(1./x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
   } else {
      y = 2/x;
      result = (std::exp(-x)/std::sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
   }
   return result;
}


/**

Probability density function of the gamma distribution.

\f[ p(x) = {1 \over \Gamma(\alpha) \theta^{\alpha}} x^{\alpha-1} e^{-x/\theta} \f]

for x>0. For detailed description see
<A HREF="http://mathworld.wolfram.com/GammaDistribution.html">
Mathworld</A>.

@ingroup PdfFunc

*/
inline double gamma_pdf(double x, double alpha, double theta, double x0 = 0) {
  // Inlined to enable clad-auto-derivation for this function.

  if ((x-x0) < 0) {
    return 0.0;
  } else if ((x-x0) == 0) {

    if (alpha == 1) {
      return 1.0/theta;
    } else {
      return 0.0;
    }

  } else if (alpha == 1) {
    return std::exp(-(x-x0)/theta)/theta;
  } else {
    return std::exp((alpha - 1) * std::log((x-x0)/theta) - (x-x0)/theta - lgam(alpha))/theta;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Computation of the normalized lower incomplete gamma function P(a,x) as defined in the
/// Handbook of Mathematical Functions by Abramowitz and Stegun, formula 6.5.1 on page 260 .
/// Its normalization is such that TMath::Gamma(a,+infinity) = 1 .
///
///  \f[
///  P(a, x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt
///  \f]
///
/// \author NvE 14-nov-1998 UU-SAP Utrecht
double Gamma(double a,double x)
{
   return igam(a, x);
}


////////////////////////////////////////////////////////////////////////////////
/// Computes the density function of Gamma distribution at point x.
///
/// \param[in] gamma   shape parameter
/// \param[in] mu      location parameter
/// \param[in] beta    scale parameter
///
/// The definition can be found in "Engineering Statistics Handbook" on site
/// http://www.itl.nist.gov/div898/handbook/eda/section3/eda366b.htm
double GammaDist(double x, double gamma, double mu, double beta)
{
   if ((x<mu) || (gamma<=0) || (beta <=0)) {
      LOG_ERROR << "tpcrs::GammaDist: illegal parameter values x = " << x << ", gamma = " << gamma << " beta = " << beta << '\n';
      return 0;
   }
   return gamma_pdf(x, gamma, beta, mu);
}

}
