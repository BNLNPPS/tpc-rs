#ifndef TPCRS_TF1F_H_
#define TPCRS_TF1F_H_

#include "TF1.h"
#include <string>

class TF1F : public TF1
{
 public:
#if ROOT_VERSION_CODE < 393216 /* = ROOT_VERSION(6,0,0) */
  TF1F() : TF1()     {fNpx       = 200;}
  TF1F(const char* name, const char* formula, double xmin = 0, double xmax = 1) :
    TF1(name, formula, xmin, xmax), fdX(-1), fStep(-1) {fNpx       = 200;}
  TF1F(const char* name, double xmin, double xmax, int npar) :
    TF1(name, xmin, xmax, npar), fdX(-1), fStep(-1) {fNpx       = 200;}
  TF1F(const char* name, void* fcn, double xmin, double xmax, int npar) :
    TF1(name, fcn, xmin, xmax, npar), fdX(-1), fStep(-1) {fNpx       = 200;}
  TF1F(const char* name, double (*fcn)(double*, double*), double xmin = 0, double xmax = 1, int npar = 0) :
    TF1(name, fcn, xmin, xmax, npar), fdX(-1), fStep(-1) {fNpx       = 200;};
#else /* ROOT 6 */
  TF1F();
  TF1F(const char* name, const char* formula, double xmin = 0, double xmax = 1);
  TF1F(const char* name, void* fcn, double xmin = 0, double xmax = 1, int npar = 0);
  TF1F(const char* name, double (*fcn)(double*, double*), double xmin = 0, double xmax = 1, int npar = 0, int ndim = 1 );
  TF1F(const char* name, double (*fcn)(const double*, const double*), double xmin = 0, double xmax = 1, int npar = 0, int ndim = 1);

  // constructor using a functor
  TF1F(const char* name, ROOT::Math::ParamFunctor f, double xmin = 0, double xmax = 1, int npar = 0, int ndim = 1);

  // Template constructors from a pointer to any C++ class of type PtrObj with a specific member function of type
  // MemFn.
  template <class PtrObj, typename MemFn>
  TF1F(const char* name, const  PtrObj &p, MemFn memFn, double xmin, double xmax, int npar, int ndim, const char* c1, const char* c2) :
    TF1(name, p, memFn, xmin, xmax, npar, c1, c2)
  {
    fNpx = 200;
  }
  // Template constructors from any  C++ callable object,  defining  the operator() (double * , double *)
  // and returning a double.
  template <typename Func>
  TF1F(const char* name, Func f, double xmin, double xmax, int npar, const char* tmp  ) :
    TF1(name, f, xmin, xmax, npar, tmp)
  {
    fNpx = 200;
  }
#endif /* ROOT 6 */
  virtual ~TF1F() {}
  virtual void Save(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
  double GetSaveL(double* xx);
  double GetSaveL(int N, double x, double* y);
 protected:
  double fXmin;
  double fXmax;
  double fdX;
  int    fStep;

};
#endif
