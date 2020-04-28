#ifndef TPCRS_DEDX_PARAMETERIZATION_H_
#define TPCRS_DEDX_PARAMETERIZATION_H_

#include "TH1.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TString.h"


enum StPidParticle {
  kPidElectron,
  kPidProton,
  kPidKaon,
  kPidPion,
  kPidMuon,
  kPidDeuteron,
  kPidTriton,
  kPidHe3,
  kPidAlpha,
  KPidParticles
};


class dEdxParameterization
{
 private:
  TString      fTag;                 //! Tag for  root file (Bichsel or PAI)
  TProfile2D*  fP;                   //! zm: The most probable value of ::log(dE/dx) versus log10(beta*gamma) and log2(dx)
  TProfile2D*  fA;                   //! mean_z: The average value of z = ::log(dE/dx) versus log10(beta*gamma) and log2(dx)
  TProfile2D*  fI70;                 //! I70: The average value after 30% truncation versus log10(beta*gamma) and log2(dx)
  TProfile2D*  fI60;                 //! I60: The average value after 40% truncation versus log10(beta*gamma) and log2(dx)
  TProfile2D*  fD;                   //! Delta_P : The most probable dE/dx versus log10(beta*gamma) and log2(dx)
  TProfile2D*  fRms;                 //! sigma_z : The RMS value of z = ::log(dE/dx) versus log10(beta*gamma) and log2(dx)
  TProfile2D*  fW;                   //! width : The RMS value of z = ::log(dE/dx) versus log10(beta*gamma) and log2(dx)
  TH3D*        fPhi;                 //! The dEdxParameterization probability versus log10(beta*gamma) and log2(dx) and z
  int        fnBins[3];            //! no. of bin for each dimension (log10(bg), log2(dx) and z = ::log(dE/dx))
  double     fbinW[3];             //! bin width
  TAxis*       fAXYZ[3];             //!
  double     fMostProbableZShift;  //!
  double     fAverageZShift;       //!
  double     fI70Shift;            //!
  double     fI60Shift;            //!
  double     fbgL10min;
  double     fbgL10max;
  double     fdxL2min;
  double     fdxL2max;
  double     fzmin;
  double     fzmax;
  TH1D*        fTrs[KPidParticles + 1][6]; //! Histograms from TpcRS dE/dx simulation for each particle type
 public :
  dEdxParameterization(const char* Tag = "p10", int keep3D = 0,
                       const double MostProbableZShift = 0,
                       const double AverageZShift      = 0,
                       const double I70Shift           = 1,
                       const double I60Shift           = 1);
  virtual ~dEdxParameterization();
  double    MostProbableZCorrection(double log10bg);
  double    I70Correction(double log10bg);
  double    GetMostProbableZ(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fMostProbableZShift + fP->Interpolate(log10bg, log2dx);
  }
  double    GetAverageZ(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fAverageZShift + MostProbableZCorrection(log10bg) + fA->Interpolate(log10bg, log2dx);
  }
  double    GetRmsZ(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fRms->Interpolate(log10bg, log2dx);
  }
  double    GetI70(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fI70Shift * fI70->Interpolate(log10bg, log2dx);
  }
  double    GetI60(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fI60Shift * fI60->Interpolate(log10bg, log2dx);
  }
  double    GetMostProbabledEdx(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fD->Interpolate(log10bg, log2dx);
  }
  double    GetdEdxWidth(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return fW->Interpolate(log10bg, log2dx);
  }
  double    GetMostProbableZM(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return MostProbableZCorrection(log10bg) + GetMostProbableZ(log10bg, log2dx);
  }
  double    GetAverageZM(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return MostProbableZCorrection(log10bg) + GetAverageZ(log10bg, log2dx);
  }
  double    GetI70M(double log10bg, double log2dx)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    return I70Correction(log10bg) * GetI70(log10bg, log2dx);
  }
  double    GetProbability(double log10bg, double log2dx, double z)
  {
    log10bg = std::max(fbgL10min, std::min(fbgL10max, log10bg));
    log2dx  = std::max(fdxL2min, std::min(fdxL2max, log2dx));
    z       = std::max(fzmin, std::min(fzmax, z));
    return fPhi->Interpolate(log10bg, log2dx, z);
  }
  void        Print();
  const char*      Tag() const {return    fTag.Data();}
  const TProfile2D*  P()   const {return     fP;}
  const TProfile2D*  A()   const {return     fA;}
  const TProfile2D*  I70() const {return   fI70;}
  const TProfile2D*  I60() const {return   fI60;}
  const TProfile2D*  D()   const {return     fD;}
  const TProfile2D*  Rms() const {return   fRms;}
  const TProfile2D*  W()   const {return     fW;}
  const TH3D*        Phi() const {return   fPhi;}
  double bgL10min() const {return fbgL10min;}
  double bgL10max() const {return fbgL10max;}
  const TH1D*       I70Trs(  int part = KPidParticles) const {return fTrs[part][0];}  // Estimation for I70 from TpcRS
  const TH1D*       I70TrsB( int part = KPidParticles) const {return fTrs[part][1];}  // Estimation for I70 - Bichsel from TpcRS
  const TH1D*       I70TrsS( int part = KPidParticles) const {return fTrs[part][2];}  // Estimation for relative sigma beta*gamma dependence for I70 from TpcRS normalized to MIP
  const TH1D*       IfitTrs( int part = KPidParticles) const {return fTrs[part][3];}  // Estimation for Ifit from TpcRS
  const TH1D*       IfitTrsB(int part = KPidParticles) const {return fTrs[part][4];}  // Estimation for Ifit - Bichsel from TpcRS
  const TH1D*       IfitTrsS(int part = KPidParticles) const {return fTrs[part][5];}  // Estimation for relative sigma beta*gamma dependence for Ifit from TpcRS normalized to MIP
  double         Get(const TH1D* hist, double log10bg) const;
  double I70Trs  (int part, double log10bg) const {return Get(fTrs[part][0], log10bg);}  // Estimation for I70 from TpcRS
  double I70TrsB (int part, double log10bg) const {return Get(fTrs[part][1], log10bg);}  // Estimation for I70 - Bichsel from TpcRS
  double I70TrsS (int part, double log10bg) const {return Get(fTrs[part][2], log10bg);}  // Estimation for relative sigma beta*gamma dependence for I70 from TpcRS normalized to MIP
  double IfitTrs (int part, double log10bg) const {return Get(fTrs[part][3], log10bg);}  // Estimation for Ifit from TpcRS
  double IfitTrsB(int part, double log10bg) const {return Get(fTrs[part][4], log10bg);}  // Estimation for Ifit - Bichsel from TpcRS
  double IfitTrsS(int part, double log10bg) const {return Get(fTrs[part][5], log10bg);}  // Estimation for relative sigma beta*gamma dependence for Ifit from TpcRS normalized to MIP
};

#endif
