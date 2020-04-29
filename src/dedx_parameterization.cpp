#include <cassert>

#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "dedx_parameterization.h"
#include "tpcrs/logger.h"

#define  PrP(A)  LOG_INFO << "\t" << (#A) << " = \t" << ( A )


dEdxParameterization::dEdxParameterization(const char* Tag, int keep3D,
    const double MostProbableZShift,
    const double AverageZShift,
    const double I70Shift,
    const double I60Shift):
  fTag (Tag), fP(0), fA(0), fI70 (0), fI60(0), fD(0),
  fRms (0), fW(0), fPhi(0),
  fMostProbableZShift(MostProbableZShift),
  fAverageZShift(AverageZShift),
  fI70Shift(I70Shift),
  fI60Shift(I60Shift),
  fbgL10min(-1), fbgL10max(4),
  fdxL2min(-0.3), fdxL2max(3),
  fzmin(-4), fzmax(6)
{
  memset (fTrs, 0, sizeof(fTrs));
  TDirectory* dir = gDirectory;
  const char*                                   rootf = "P10T.root";

  if (fTag.Contains("pai", TString::kIgnoreCase)) rootf = "PaiT.root";

  if (fTag.Contains("p10", TString::kIgnoreCase)) rootf = "P10T.root";

  if (fTag.Contains("bich", TString::kIgnoreCase)) rootf = "BichselT.root";

  static const char* path  = ".:./StarDb/dEdxModel:./StarDb/global/dEdx:./StRoot/StBichsel:$STAR/StarDb/dEdxModel:$STAR/StarDb/global/dEdx:$STAR/StRoot/StBichsel";
  char* file = gSystem->Which(path, rootf, kReadPermission);

  assert(file && "File not found");

  TFile*       pFile = new TFile(file);
  delete [] file;
  assert(pFile);
  fP   = (TProfile2D*) pFile->Get("bichP");   assert(fP);   fP->SetDirectory(0);
  fA   = (TProfile2D*) pFile->Get("bichA");   assert(fA);   fA->SetDirectory(0);
  fI70 = (TProfile2D*) pFile->Get("bichI70"); assert(fI70); fI70->SetDirectory(0);
  fI60 = (TProfile2D*) pFile->Get("bichI60"); assert(fI60); fI60->SetDirectory(0);
  fD   = (TProfile2D*) pFile->Get("bichD");   assert(fD);   fD->SetDirectory(0);
  fRms = (TProfile2D*) pFile->Get("bichRms"); assert(fRms); fRms->SetDirectory(0);
  fW   = (TProfile2D*) pFile->Get("bichW");   assert(fW);   fW->SetDirectory(0);
  fPhi = (TH3D*) pFile->Get("bichPhi"); assert(fPhi); fPhi->SetDirectory(0);
  fbgL10min = fPhi->GetXaxis()->GetBinCenter(1) + 1e-7;
  fbgL10max = fPhi->GetXaxis()->GetBinCenter(fPhi->GetXaxis()->GetNbins()) - 1e-7;
  fdxL2min  = fPhi->GetYaxis()->GetBinCenter(1) + 1e-7;
  fdxL2max  = fPhi->GetYaxis()->GetBinCenter(fPhi->GetYaxis()->GetNbins()) - 1e-7;
  fzmin  = fPhi->GetZaxis()->GetBinCenter(1) + 1e-7;
  fzmax  = fPhi->GetZaxis()->GetBinCenter(fPhi->GetZaxis()->GetNbins()) - 1e-7;

  if (dir) dir->cd();

  for (int i = 0; i < 3; i++) {
    if (i == 0) fAXYZ[i] = fPhi->GetXaxis();

    if (i == 1) fAXYZ[i] = fPhi->GetYaxis();

    if (i == 2) fAXYZ[i] = fPhi->GetZaxis();

    fnBins[i] = fAXYZ[i]->GetNbins();
    fbinW[i]  = fAXYZ[i]->GetBinWidth(1);
    PrP(i); PrP(fnBins[i]); PrP(fbinW[i]); LOG_INFO << '\n';
    assert(fnBins[i] != 1);
  }

  //  if (! keep3D) SafeDelete(fPhi);
  // set normalization factor to 2.3976 keV/cm at beta*gamma = 4;
  static const double dEdxMIP = 2.39761562607903311; // [keV/cm]
  static const double MIPBetaGamma10 = std::log10(4.);
  //  fMostProbableZShift = std::log(dEdxMIP) - Interpolation(fP,MIPBetaGamma10,1,0);
  //  fAverageZShift      = std::log(dEdxMIP) - Interpolation(fA,MIPBetaGamma10,1,0);
  fI70Shift           *= dEdxMIP / GetI70(MIPBetaGamma10, 1);
  fI60Shift           *= dEdxMIP / GetI60(MIPBetaGamma10, 1);
  fMostProbableZShift  = std::log(fI70Shift);
  fAverageZShift       = fMostProbableZShift;
  const char* Names[KPidParticles + 1] = {"e", "proton", "kaon", "pi", "mu", "deuteron", "triton", "He3", "alpha", "all"};

  for (int i = 0; i <= KPidParticles; i++) {
    TString name(Names[i]);
    const char* type[6] = {"70p", "70", "70S", "zp", "z", "zS"};

    for (int j = 0; j < 6; j++) {
      fTrs[i][j] = (TH1D*) pFile->Get(name + type[j]);

      if (fTrs[i][j])  fTrs[i][j]->SetDirectory(0);
    }
  }

  delete pFile;
}


dEdxParameterization::~dEdxParameterization()
{
  SafeDelete(fP);
  SafeDelete(fA);
  SafeDelete(fI70);
  SafeDelete(fI60);
  SafeDelete(fD);
  SafeDelete(fRms);
  SafeDelete(fW);
  SafeDelete(fPhi);

  for (int i = 0; i <= KPidParticles; i++)
    for (int j = 0; j < 6; j++) {SafeDelete(fTrs[i][j]);}
}


void dEdxParameterization::Print()
{
  PrP(fTag); LOG_INFO << '\n';
  PrP(fP); if (fP) PrP(fP->GetTitle()); LOG_INFO << '\n';
  PrP(fA); if (fA) PrP(fA->GetTitle()); LOG_INFO << '\n';
  PrP(fI70); if (fI70) PrP(fI70->GetTitle()); LOG_INFO << '\n';
  PrP(fI60); if (fI60) PrP(fI60->GetTitle()); LOG_INFO << '\n';
  PrP(fD); if (fD) PrP(fD->GetTitle()); LOG_INFO << '\n';
  PrP(fRms); if (fRms) PrP(fRms->GetTitle()); LOG_INFO << '\n';
  PrP(fW); if (fW) PrP(fW->GetTitle()); LOG_INFO << '\n';
  PrP(fPhi); if (fPhi) PrP(fPhi->GetTitle()); LOG_INFO << '\n';
  PrP(fMostProbableZShift); LOG_INFO << '\n';
  PrP(fAverageZShift); LOG_INFO << '\n';
  PrP(fI70Shift); LOG_INFO << '\n';
  PrP(fI60Shift); LOG_INFO << '\n';
}


double dEdxParameterization::MostProbableZCorrection(double log10bg)
{
  static const double pars[2] = {-3.68846e-03, 4.72944e+00}; // FitHzAllHist012P05id  FitH + Prof 050905
  return pars[0] * std::exp(-pars[1] * log10bg);
}

double dEdxParameterization::I70Correction(double log10bg)
{
  static const double pars[2] = {-1.65714e-02, 3.27271e+00}; //  FitH70AllHist012P05id FitH + Prof 050905
  return std::exp(pars[0] * std::exp(-pars[1] * log10bg));
}


double dEdxParameterization::Get(const TH1D* hist, double log10bg) const
{
  static TH1D* hsave = 0;
  static double xmin = -100, xmax = 100;

  if (hist != hsave) {
    hsave = (TH1D*) hist;
    TAxis* x = hsave->GetXaxis();
    int f = x->GetFirst();
    int l = x->GetLast();

    xmin = x->GetBinUpEdge(f);
    xmax = x->GetBinLowEdge(l);
  }

  if (log10bg < xmin) log10bg = xmin;

  if (log10bg > xmax) log10bg = xmax;

  return hsave->Interpolate(log10bg);
}
