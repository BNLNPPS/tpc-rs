#pragma once

namespace TCL {

#define TCL_MXMAD(n_,a,b,c,i,j,k)                       \
     /* Local variables */                                \
     int l, m, n, ia, ic, ib, ja, jb, iia, iib, ioa, iob; \
                                                          \
     /* Parameter adjustments */                          \
     --a;  --b;  --c;                                     \
     /* Function Body */                                  \
 /*                      MXMAD MXMAD1 MXMAD2 MXMAD3 MXMPY MXMPY1 MXMPY2 MXMPY3 MXMUB MXMUB1 MXMUB2 MXMUB3 */ \
 /*  const int iandj1[] = {21,   22,    23,    24,   11,    12,    13,    14,    31,   32,   33,    34 }; */ \
     const int iandj1[] = {2,    2 ,    2 ,    2 ,   1 ,    1 ,    1 ,    1 ,    3 ,   3 ,   3 ,    3  }; \
     const int iandj2[] = { 1,    2,     3,     4,    1,     2,     3,     4,     1,    2,    3,     4 }; \
     int n1 = iandj1[n_];                                  \
     int n2 = iandj2[n_];                                  \
     if (i == 0 || k == 0) return 0;                       \
                                                           \
     switch (n2) {                                         \
       case 1: iia = 1; ioa = j; iib = k; iob = 1; break;  \
       case 2: iia = 1; ioa = j; iib = 1; iob = j; break;  \
       case 3: iia = i; ioa = 1; iib = k; iob = 1; break;  \
       case 4: iia = i; ioa = 1; iib = 1; iob = j; break;  \
       default: iia = ioa = iib = iob = 0; assert(iob);    \
     };                                                    \
                                                           \
     ia = 1; ic = 1;                                       \
     for (l = 1; l <= i; ++l) {                            \
             ib = 1;                                           \
             for (m = 1; m <= k; ++m,++ic) {                   \
               switch (n1) {                                   \
                       case 1:  c[ic] = 0.;      break;            \
                       case 3:  c[ic] = -c[ic];  break;            \
               };                                              \
               if (j == 0) continue;                           \
               ja = ia; jb = ib;                               \
           double cic = c[ic];                             \
               for (n = 1; n <= j; ++n, ja+=iia, jb+=iib)      \
                        cic += a[ja] * b[jb];                      \
           c[ic] = cic;                                    \
               ib += iob;                                      \
             }                                                 \
             ia += ioa;                                        \
     }

double *mxmad_0_(int n_, const double *a, const double *b, double *c, int i, int j, int k)
{
  TCL_MXMAD(n_,a,b,c,i,j,k)
  return c;
}

// TCL::mxmpy(A.GetArray(), B.GetArray(), fArray, NI, NJ, NK);
// NI = A.GetNrows(); fNcols = NI;
// NJ = A.GetNcols();
// NK = B.GetNcols(); fNrows = NK;
inline double *mxmpy(const double *a, const double *b, double *c, int i, int j, int k)
{
  return mxmad_0_(4, a, b, c, i, j, k);
}

}
