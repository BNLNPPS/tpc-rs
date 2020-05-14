#ifndef TPCRS_BICHSEL_H_
#define TPCRS_BICHSEL_H_

#include "TString.h"

#include "dedx_parameterization.h"


class Bichsel
{
 public:
  enum EParTypes {kP10, kBichsel, kPAI, kTotal};
 private:
  static TString        m_Tags[kTotal];
  int                 m_Type;
  TString               m_Tag;
  static dEdxParameterization* m_dEdxParameterizations[kTotal]; //!
  dEdxParameterization* m_dEdxParameterization; //!
  static Bichsel*       fgBichsel; //! last instance
 public:
  Bichsel(const char* tag = "P10", int keep3D = 0);
  ~Bichsel() {fgBichsel = 0;};
  static Bichsel* Instance(const char* tag = "P10", int keep3D = 0);
  static void Clean();
  double    GetMostProbableZ(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetMostProbableZ(log10bg, log2dx);
  }
  double    GetMostProbableZM(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetMostProbableZM(log10bg, log2dx);
  }
  double    GetAverageZ(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetAverageZ(log10bg, log2dx);
  }
  double    GetAverageZM(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetAverageZM(log10bg, log2dx);
  }
  double    GetRmsZ(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetRmsZ(log10bg, log2dx);
  }
  double    GetI70(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetI70(log10bg, log2dx);
  }
  double    GetI70M(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetI70M(log10bg, log2dx);
  }
  double    GetI60(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetI60(log10bg, log2dx);
  }
  double    GetMostProbabledEdx(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetMostProbabledEdx(log10bg, log2dx);
  }
  double    GetdEdxWidth(double log10bg, double log2dx = 1.)
  {
    return m_dEdxParameterization->GetdEdxWidth(log10bg, log2dx);
  }
  double    GetProbability(double log10bg, double log2dx, double z)
  {
    return m_dEdxParameterization->GetProbability(log10bg, log2dx, z);
  }
  const dEdxParameterization* Parameterization() const {return m_dEdxParameterization;}
  void Print();
  const char*      Tag() const {return    m_dEdxParameterization->Tag();}
  const TProfile2D*  P()   const {return     m_dEdxParameterization->P();}
  const TProfile2D*  A()   const {return     m_dEdxParameterization->A();}
  const TProfile2D*  I70() const {return   m_dEdxParameterization->I70();}
  const TProfile2D*  I60() const {return   m_dEdxParameterization->I60();}
  const TProfile2D*  D()   const {return     m_dEdxParameterization->D();}
  const TProfile2D*  Rms() const {return   m_dEdxParameterization->Rms();}
  const TProfile2D*  W()   const {return     m_dEdxParameterization->W();}
  const TH3D*        Phi() const {return   m_dEdxParameterization->Phi();}
  const TH1D*       I70Trs  (int part = KPidParticles) const {return m_dEdxParameterization->I70Trs( part);}  // Estimation for I70 from TpcRS
  const TH1D*       I70TrsB (int part = KPidParticles) const {return m_dEdxParameterization->I70TrsB(part);}  // Estimation for I70 - Bichsel from TpcRS
  const TH1D*       I70TrsS (int part = KPidParticles) const {return m_dEdxParameterization->I70TrsS(part);}  // Estimation for relative sigma bg dependence for I70 from TpcRS normalized to MIP
  const TH1D*       IfitTrs (int part = KPidParticles) const {return m_dEdxParameterization->IfitTrs(part);}  // Estimation for Ifit from TpcRS
  const TH1D*       IfitTrsB(int part = KPidParticles) const {return m_dEdxParameterization->IfitTrsB(part);} // Estimation for Ifit - Bichsel from TpcRS
  const TH1D*       IfitTrsS(int part = KPidParticles) const {return m_dEdxParameterization->IfitTrsS(part);} // Estimation for relative sigma bg dependence for Ifit from TpcRS normalized to MIP
  double I70Trs  (int part, double log10bg) const {return m_dEdxParameterization->Get(I70Trs  (part), log10bg);}  // Estimation for I70 from TpcRS
  double I70TrsB (int part, double log10bg) const {return m_dEdxParameterization->Get(I70TrsB (part), log10bg);}  // Estimation for I70 - Bichsel from TpcRS
  double I70TrsS (int part, double log10bg) const {return m_dEdxParameterization->Get(I70TrsS (part), log10bg);}  // Estimation for relative sigma beta*gamma dependence for I70 from TpcRS normalized to MIP
  double IfitTrs (int part, double log10bg) const {return m_dEdxParameterization->Get(IfitTrs (part), log10bg);}  // Estimation for Ifit from TpcRS
  double IfitTrsB(int part, double log10bg) const {return m_dEdxParameterization->Get(IfitTrsB(part), log10bg);}  // Estimation for Ifit - Bichsel from TpcRS
  double IfitTrsS(int part, double log10bg) const {return m_dEdxParameterization->Get(IfitTrsS(part), log10bg);}  // Estimation for relative sigma beta*gamma dependence for Ifit from TpcRS normalized to MIP
};

#endif
