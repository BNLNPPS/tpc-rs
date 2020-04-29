#include "TF1F.h"
#include "math_funcs.h"


void TF1F::Save(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  fXmin = xmin; fXmax = xmax; fStep = 20, fdX = 1. / fStep; fNpx = tpcrs::irint((fXmax - fXmin) / fdX);
  TF1::Save(xmin, xmax, ymin, ymax, zmin, zmax);
}


double TF1F::GetSaveL(double* xx)
{
  // Get value corresponding to X in array of fSave values
  if (xx[0] < fXmin || xx[0]  > fXmax || fdX <= 0) return 0.;

  int bin     = tpcrs::irint((xx[0] - fXmin) / fdX);
  return fSave[bin];
}


double TF1F::GetSaveL(int N, double x, double* y)
{
  // Get values y[N] corresponding to x+i, i = [0, ..., N-1];
  //  memset(y, 0, N*sizeof(double));
  int bin     = tpcrs::irint((x - fXmin) / fdX);
  int i1 = 0;

  while (bin < 0) {i1++; bin += fStep;}

  for (int i = i1; i < N && bin < GetNpx() - 3; i++, bin += fStep) {
    y[i] = fSave[bin];
  }

  return y[0];
}


double TF1F::GetSaveL(int N, double* x, double* y)
{
  // Get values y[N] corresponding to x[N] in array of fSave values
  memset(y, 0, N * sizeof(double));

  if (GetNpx() <= 0) return 0.;

  for (int i = 0; i < N; i++) {
    if (x[i] > fXmin) {
      int bin     = int((x[i] - fXmin) / fdX);

      if (bin < GetNpx() - 3) y[i] = fSave[bin];
    }
  }

  return y[0];
}
#if ROOT_VERSION_CODE >= 393216 /* = ROOT_VERSION(6,0,0) */
TF1F::TF1F(): TF1() {fNpx = 200;}
TF1F::TF1F(const char* name, const char* formula, double xmin, double xmax)
  : TF1(name, formula, xmax, xmin) {fNpx = 200;}
TF1F::TF1F(const char* name, double (*fcn)(double*, double*), double xmin, double xmax, int npar, int ndim)
  : TF1(name, fcn, xmin, xmax, npar, ndim) {fNpx = 200;}
#endif /* ROOT 6 */
