/**
 *  The maker's algorithms and formulae based on
 *  http://www.inst.bnl.gov/programs/gasnobledet/publications/Mathieson's_Book.pdf,
 *  and  Photo Absorption Model
 *  "Ionization energy loss in very thin absorbers.", V.M. Grishin, V.K. Ermilova, S.K. Kotelnikov Nucl.Instrum.Meth.A309:476-484,1991
 *  "A method to improve tracking and particle identification in TPC's and silicon detectors.", Hans Bichsel, Nucl.Instrum.Meth.A562:154-197,2006
 *  HEED: "Modeling of ionization produced by fast charged particles in gases", I.B. Smirnov, Nucl.Instrum.Meth.A55(2005) 747-493.
 *
 *  \Author Y.Fisyak, fisyak@bnl.gov
 */

#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

#include "StTpcRSMaker/StTpcRSMaker.h"
// SCL
#include "StarClassLibrary/StThreeVector.hh"
#include "StarClassLibrary/StPhysicalHelix.hh"
// ROOT
#include "TTableSorter.h"
#include "TRandom.h"
#include "TFile.h"
#include "TBenchmark.h"
#include "TProfile2D.h"
#include "StTpcRSMaker/tcl.h"
#include "Math/SpecFuncMathMore.h"

#include "StDbUtilities/StCoordinates.hh"
#include "StDbUtilities/StTpcCoordinateTransform.hh"
#include "StDbUtilities/StMagUtilities.h"
#include "StDetectorDbMaker/St_tpcAltroParamsC.h"
#include "StDetectorDbMaker/St_asic_thresholdsC.h"
#include "StDetectorDbMaker/St_tss_tssparC.h"
#include "StDetectorDbMaker/St_tpcPadGainT0BC.h"
#include "StDetectorDbMaker/St_TpcResponseSimulatorC.h"
#include "StDetectorDbMaker/St_tpcAnodeHVavgC.h"
#include "StDetectorDbMaker/StDetectorDbTpcRDOMasks.h"
#include "StDetectorDbMaker/St_tpcPadConfigC.h"
#include "StDetectorDbMaker/St_tpcGainCorrectionC.h"
#include "StDetectorDbMaker/St_TpcAvgCurrentC.h"
#include "StDetectorDbMaker/St_TpcAvgPowerSupplyC.h"
#include "StDetectorDbMaker/St_trigDetSumsC.h"
#include "StarClassLibrary/StParticleTable.hh"
#include "StarClassLibrary/StParticleDefinition.hh"
#include "StTpcRSMaker/Altro.h"
#include "StBichsel/Bichsel.h"
#include "StdEdxY2Maker/StTpcdEdxCorrection.h"

#include "tpcrs/configurator.h"
#include "tpcrs/logger.h"
#include "tpcrs/math.h"

struct HitPoint_t {
  int indx;
  int TrackId;
  double s; // track length to current point
  double sMin, sMax;
  g2t_tpc_hit_st* tpc_hitC;
  StGlobalCoordinate   xyzG;
  StTpcLocalSectorCoordinate coorLS;
  StTpcLocalSectorDirection dirLS, BLS;
  StTpcPadCoordinate Pad;
};

#define __STOPPED_ELECTRONS__
#define __DEBUG__
#if defined(__DEBUG__)
#define PrPP(A,B) if (Debug()%10 > 2) {LOG_INFO << "StTpcRSMaker::" << (#A) << "\t" << (#B) << " = \t" << (B) << '\n';}
#else
#define PrPP(A,B)
#endif
static bool ClusterProfile = false;
#define Laserino 170
#define Chasrino 171

//                                    Inner        Outer
static       double t0IO[2]   = {1.20868e-9, 1.43615e-9}; // recalculated in InducedCharge
static const double tauC[2]   = {999.655e-9, 919.183e-9};
TF1F*     StTpcRSMaker::fgTimeShape3[2]    = {0, 0};
TF1F*     StTpcRSMaker::fgTimeShape0[2]    = {0, 0};

static const int nx[2] = {200, 500};
static const double xmin[2] =  {-10., -6};
static const double xmax[2] =  { 10., 44};
static const int nz = 42;
static const double zmin = -210;
static const double zmax = -zmin;
//                     io pt
static TProfile2D* hist[5][3] = {0};
static const int nChecks = 21;
static TH1*  checkList[2][21] = {0};
static TString TpcMedium("TPCE_SENSITIVE_GAS");


StTpcRSMaker::StTpcRSMaker(double eCutOff, const char* name):
  mLaserScale(1),
  minSignal(1e-4),
  ElectronRange(0.0055), // Electron Range(.055mm)
  ElectronRangeEnergy(3000), // eV
  ElectronRangePower(1.78), // sigma =  ElectronRange*(eEnery/ElectronRangeEnergy)**ElectronRangePower
  max_sectors(24),
  max_pads(182),
  mCutEle(eCutOff)
{
  m_Mode = 0;
  //  SETBIT(m_Mode,kHEED);
  SETBIT(m_Mode, kBICHSEL); // Default is Bichsel
  SETBIT(m_Mode, kdEdxCorr);
  SETBIT(m_Mode, kDistortion);
  InitRun(0);
}


StTpcRSMaker::~StTpcRSMaker()
{
  SafeDelete(mAltro);
  Finish();
}


void StTpcRSMaker::Finish()
{
  if (m_SignalSum) {free(m_SignalSum); m_SignalSum = 0;}

  SafeDelete(mdNdx);
  SafeDelete(mdNdxL10);
  SafeDelete(mdNdEL10);

  for (int io = 0; io < 2; io++) {// Inner/Outer
    for (int sec = 0; sec < max_sectors; sec++) {
      if (mShaperResponses[io][sec] && !mShaperResponses[io][sec]->TestBit(TObject::kNotDeleted)) {SafeDelete(mShaperResponses[io][sec]);}

      SafeDelete(mChargeFraction[io][sec]);
      SafeDelete(mPadResponseFunction[io][sec]);
    }

    SafeDelete(mPolya[io]);
  }

  if (m_TpcdEdxCorrection && m_TpcdEdxCorrection->TestBit(kCanDelete)) delete m_TpcdEdxCorrection;

  m_TpcdEdxCorrection = 0;
}


void StTpcRSMaker::InitRun(int runnumber)
{
  if (!gStTpcDb) {
    LOG_ERROR << "Database Missing! Can't initialize TpcRS\n";
    return;
  }

  if (TESTBIT(m_Mode, kBICHSEL)) {
    LOG_INFO << "StTpcRSMaker:: use H.Bichsel model for dE/dx simulation\n";

    TFile inner(tpcrs::Configurator::Locate("dNdE_Bichsel.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(tpcrs::Configurator::Locate("dNdx_Bichsel.root").c_str());
    mdNdx = (TH1D*) outer.Get("dNdx"); assert(mdNdx);
    mdNdx->SetDirectory(0);
  }
  else if (TESTBIT(m_Mode, kHEED)) {
    LOG_INFO << "StTpcRSMaker:: use Heed model for dE/dx simulation\n";

    TFile inner(tpcrs::Configurator::Locate("dNdx_Heed.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(tpcrs::Configurator::Locate("dNdx_Heed.root").c_str());
    mdNdxL10 = (TH1D*) outer.Get("dNdxL10"); assert(mdNdxL10);
    mdNdxL10->SetDirectory(0);
  }
  else {LOG_INFO << "StTpcRSMaker:: use GEANT321 model for dE/dx simulation\n";}

  // Distortions
  if (TESTBIT(m_Mode, kdEdxCorr)) {
    LOG_INFO << "StTpcRSMaker:: use Tpc dE/dx correction from calibaration\n";
    int Mask = -1; // 22 bits
    CLRBIT(Mask, StTpcdEdxCorrection::kAdcCorrection);
    CLRBIT(Mask, StTpcdEdxCorrection::kAdcCorrectionMDF);
    CLRBIT(Mask, StTpcdEdxCorrection::kdXCorrection);
    //    CLRBIT(Mask,StTpcdEdxCorrection::kEdge);
    //    CLRBIT(Mask,StTpcdEdxCorrection::kTanL);
    m_TpcdEdxCorrection = new StTpcdEdxCorrection(Mask, Debug());
  }

  if (TESTBIT(m_Mode, kDistortion)) {
    LOG_INFO << "StTpcRSMaker:: use Tpc distortion correction\n";
  }

  double samplingFrequency     = 1.e6 * gStTpcDb->Electronics()->samplingFrequency(); // Hz
  double TimeBinWidth          = 1. / samplingFrequency;
  /*
  select firstInnerSectorAnodeWire,lastInnerSectorAnodeWire,numInnerSectorAnodeWires,firstOuterSectorAnodeWire,lastOuterSectorAnodeWire,numOuterSectorAnodeWires from  Geometry_tpc.tpcWirePlanes;
  +---------------------------+--------------------------+--------------------------+---------------------------+--------------------------+--------------------------+
  | firstInnerSectorAnodeWire | lastInnerSectorAnodeWire | numInnerSectorAnodeWires | firstOuterSectorAnodeWire | lastOuterSectorAnodeWire | numOuterSectorAnodeWires |
  +---------------------------+--------------------------+--------------------------+---------------------------+--------------------------+--------------------------+
  |             53.2000000000 |           120.8000000000 |                      170 |            122.7950000000 |           191.1950000000 |                      172 |
  +---------------------------+--------------------------+--------------------------+---------------------------+--------------------------+--------------------------+
   */
  numberOfInnerSectorAnodeWires  = gStTpcDb->WirePlaneGeometry()->numberOfInnerSectorAnodeWires ();
  firstInnerSectorAnodeWire      = gStTpcDb->WirePlaneGeometry()->firstInnerSectorAnodeWire();
  lastInnerSectorAnodeWire       = gStTpcDb->WirePlaneGeometry()->lastInnerSectorAnodeWire ();
  numberOfOuterSectorAnodeWires  = gStTpcDb->WirePlaneGeometry()->numberOfOuterSectorAnodeWires ();
  firstOuterSectorAnodeWire      = gStTpcDb->WirePlaneGeometry()->firstOuterSectorAnodeWire();
  lastOuterSectorAnodeWire       = gStTpcDb->WirePlaneGeometry()->lastOuterSectorAnodeWire ();
  anodeWirePitch                 = gStTpcDb->WirePlaneGeometry()->anodeWirePitch           ();
  anodeWireRadius                = gStTpcDb->WirePlaneGeometry()->anodeWireRadius();
  float BFieldG[3];
  float xyz[3] = {0, 0, 0};
  StarMagField::Instance()->BField(xyz, BFieldG);
  // Shapers
  double timeBinMin = -0.5;
  double timeBinMax = 44.5;
  const char* Names[2] = {"I", "O"};
  double CathodeAnodeGap[2] = {0.2, 0.4};

  for (int sector = 1; sector <= 24; sector++) {
    innerSectorAnodeVoltage[sector - 1] = outerSectorAnodeVoltage[sector - 1] = 0;
    int nAliveInner = 0;
    int nAliveOuter = 0;

    for (int row = 1; row <= St_tpcPadConfigC::instance()->numberOfRows(sector); row++) {
      if (St_tpcPadConfigC::instance()->IsRowInner(sector, row)) {
        nAliveInner++;
        innerSectorAnodeVoltage[sector - 1] += St_tpcAnodeHVavgC::instance()->voltagePadrow(sector, row);
      }
      else {
        nAliveOuter++;
        outerSectorAnodeVoltage[sector - 1] += St_tpcAnodeHVavgC::instance()->voltagePadrow(sector, row);
      }
    }

    if (! nAliveInner && ! nAliveOuter) {
      LOG_INFO << "Illegal date/time. Tpc sector " << sector << " Anode Voltage is not set to run condition: AliveInner: " << nAliveInner
               << "\tAliveOuter: " << nAliveOuter
               << "\tStop the run\n";
      assert(nAliveInner || nAliveOuter);
    }
    else {
      if (nAliveInner > 1) innerSectorAnodeVoltage[sector - 1] /= nAliveInner;

      if (nAliveOuter > 1) outerSectorAnodeVoltage[sector - 1] /= nAliveOuter;
    }

    for (int io = 0; io < 2; io++) {// In/Out
      if (io == 0) {
        if (sector > 1 && TMath::Abs(innerSectorAnodeVoltage[sector - 1] - innerSectorAnodeVoltage[sector - 2]) < 1) {
          InnerAlphaVariation[sector - 1] = InnerAlphaVariation[sector - 2];
        }
        else {
          LOG_INFO << "Inner Sector " << sector << " ======================\n";
          InnerAlphaVariation[sector - 1] = InducedCharge(anodeWirePitch,
                                            CathodeAnodeGap[io],
                                            anodeWireRadius,
                                            innerSectorAnodeVoltage[sector - 1], t0IO[io]);
        }
      }
      else {
        if (sector > 1 && TMath::Abs(outerSectorAnodeVoltage[sector - 1] - outerSectorAnodeVoltage[sector - 2]) < 1) {
          OuterAlphaVariation[sector - 1] = OuterAlphaVariation[sector - 2];
        }
        else {
          LOG_INFO << "Outer Sector " << sector << " ======================\n";
          OuterAlphaVariation[sector - 1] = InducedCharge(anodeWirePitch,
                                            CathodeAnodeGap[io],
                                            anodeWireRadius,
                                            outerSectorAnodeVoltage[sector - 1], t0IO[io]);
        }
      }
    }
  }

  for (int io = 0; io < 2; io++) {// In/Out
    //  mPolya = new TF1F("Polya;x = G/G_0;signal","sqrt(x)/exp(1.5*x)",0,10); // original Polya
    //  mPolya = new TF1F("Polya;x = G/G_0;signal","pow(x,0.38)*exp(-1.38*x)",0,10); //  Valeri Cherniatin
    //   mPoly = new TH1D("Poly","polyaAvalanche",100,0,10);
    double gamma;

    if (!io ) gamma = St_TpcResponseSimulatorC::instance()->PolyaInner();
    else      gamma = St_TpcResponseSimulatorC::instance()->PolyaOuter();

    if (gamma <= 0) gamma = 1.38;

    mPolya[io] = new TF1F(io == 0 ? "PolyaInner;x = G/G_0;signal" : "PolyaOuter;x = G/G_0;signal", polya, 0, 10, 3);
    mPolya[io]->SetParameters(gamma, 0., 1. / gamma);
    double params3[7] = {t0IO[io],
                           St_TpcResponseSimulatorC::instance()->tauF(),
                           St_TpcResponseSimulatorC::instance()->tauP(),
                           St_TpcResponseSimulatorC::instance()->tauIntegration(),
                           TimeBinWidth,     0, (double ) io
                          };
    double params0[5] = {t0IO[io], St_TpcResponseSimulatorC::instance()->tauX()[io], TimeBinWidth,     0, (double ) io};

    if (! fgTimeShape3[io]) {// old electronics, intergation + shaper alltogether
      fgTimeShape3[io] = new TF1F(Form("TimeShape3%s", Names[io]),
                                  shapeEI3, timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth, 7);
      fgTimeShape3[io]->SetParNames("t0", "tauF", "tauP", "tauI", "width", "tauC", "io");
      fgTimeShape3[io]->SetParameters(params3);
      params3[5] = fgTimeShape3[io]->Integral(timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth);
      fgTimeShape3[io]->SetTitle(fgTimeShape3[io]->GetName());
      fgTimeShape3[io]->GetXaxis()->SetTitle("time (secs)");
      fgTimeShape3[io]->GetYaxis()->SetTitle("signal");
    }

    if (! fgTimeShape0[io]) {// new electronics only integration
      fgTimeShape0[io] = new TF1F(Form("TimeShape%s", Names[io]),
                                  shapeEI, timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth, 5);
      fgTimeShape0[io]->SetParNames("t0", "tauI", "width", "tauC", "io");
      params0[3] = St_TpcResponseSimulatorC::instance()->tauC()[io];
      fgTimeShape0[io]->SetParameters(params0);
      params0[3] = fgTimeShape0[io]->Integral(0, timeBinMax * TimeBinWidth);
      fgTimeShape0[io]->SetTitle(fgTimeShape0[io]->GetName());
      fgTimeShape0[io]->GetXaxis()->SetTitle("time (secs)");
      fgTimeShape0[io]->GetYaxis()->SetTitle("signal");
    }

    for (int sector = 1; sector <= max_sectors; sector++) {
      //                             w       h         s      a       l   i
      //  double paramsI[6] = {0.2850, 0.2000,  0.4000, 0.0010, 1.1500, 0};
      //  double paramsO[6] = {0.6200, 0.4000,  0.4000, 0.0010, 1.1500, 0};
      double xmaxP =  4.5;//4.5*St_tpcPadConfigC::instance()->innerSectorPadWidth(sector);// 4.5
      double xminP = -xmaxP;
      double params[6];

      if (! io) {
        params[0] = St_tpcPadConfigC::instance()->innerSectorPadWidth(sector);               // w = width of pad
        params[1] = gStTpcDb->WirePlaneGeometry()->innerSectorAnodeWirePadPlaneSeparation(); // h = Anode-Cathode gap
        params[2] = gStTpcDb->WirePlaneGeometry()->anodeWirePitch();                         // s = wire spacing
        params[3] = St_TpcResponseSimulatorC::instance()->K3IP();
        params[4] = 0;
        params[5] = St_tpcPadConfigC::instance()->innerSectorPadPitch(sector);
      }
      else {
        params[0] = St_tpcPadConfigC::instance()->outerSectorPadWidth(sector);              // w = width of pad
        params[1] = gStTpcDb->WirePlaneGeometry()->outerSectorAnodeWirePadPlaneSeparation();// h = Anode-Cathode gap
        params[2] = gStTpcDb->WirePlaneGeometry()->anodeWirePitch();                        // s = wire spacing
        params[3] = St_TpcResponseSimulatorC::instance()->K3OP();
        params[4] = 0;
        params[5] = St_tpcPadConfigC::instance()->outerSectorPadPitch(sector);
      }

      if ( !mPadResponseFunction[io][sector - 1] ) {
        mPadResponseFunction[io][sector - 1] = new TF1F(io == 0 ? "PadResponseFunctionInner" : "PadResponseFunctionOuter", StTpcRSMaker::PadResponseFunc, xminP, xmaxP, 6);
        mPadResponseFunction[io][sector - 1]->SetParameters(params);
        mPadResponseFunction[io][sector - 1]->SetParNames("PadWidth", "Anode-Cathode gap", "wire spacing", "K3OP", "CrossTalk", "PadPitch");
        mPadResponseFunction[io][sector - 1]->SetTitle(mPadResponseFunction[io][sector - 1]->GetName());
        mPadResponseFunction[io][sector - 1]->GetXaxis()->SetTitle("pads");
        mPadResponseFunction[io][sector - 1]->GetYaxis()->SetTitle("Signal");
        // Cut tails
        double x = 2.5;//cm   xmaxP;
        mPadResponseFunction[io][sector - 1]->SetRange(-x, x);
        mPadResponseFunction[io][sector - 1]->Save(xminP, xmaxP, 0, 0, 0, 0);
      }

      if ( !mChargeFraction[io][sector - 1] ) {
        xmaxP = 2.5;//5*St_tpcPadConfigC::instance()->innerSectorPadLength(sector); // 1.42
        xminP = -xmaxP;
        mChargeFraction[io][sector - 1] = new TF1F(io == 0 ? "ChargeFractionInner" : "ChargeFractionOuter", StTpcRSMaker::PadResponseFunc, xminP, xmaxP, 6);

        if (!io) {
          params[0] = St_tpcPadConfigC::instance()->innerSectorPadLength(sector);
          params[3] = St_TpcResponseSimulatorC::instance()->K3IR();
          params[4] = 0;
          params[5] = 1.;
        }
        else {
          params[0] = St_tpcPadConfigC::instance()->outerSectorPadLength(sector);
          params[3] = St_TpcResponseSimulatorC::instance()->K3OR();
          params[4] = 0;
          params[5] = 1.;
        }

        mChargeFraction[io][sector - 1]->SetParameters(params);
        mChargeFraction[io][sector - 1]->SetParNames("PadLength", "Anode-Cathode gap", "wire spacing", "K3IR", "CrossTalk", "RowPitch");
        mChargeFraction[io][sector - 1]->SetTitle(mChargeFraction[io][sector - 1]->GetName());
        mChargeFraction[io][sector - 1]->GetXaxis()->SetTitle("Distance (cm)");
        mChargeFraction[io][sector - 1]->GetYaxis()->SetTitle("Signal");
        // Cut tails
        double x = xmaxP;
        double ymax = mChargeFraction[io][sector - 1]->Eval(0);

        for (; x > 1.5; x -= 0.05) {
          double r = mChargeFraction[io][sector - 1]->Eval(x) / ymax;

          if (r > 1e-2) break;
        }

        mChargeFraction[io][sector - 1]->SetRange(-x, x);
        mChargeFraction[io][sector - 1]->Save(xminP, xmaxP, 0, 0, 0, 0);
      }

      //  TF1F *func = new TF1F("funcP","x*sqrt(x)/exp(2.5*x)",0,10);
      // see http://www4.rcf.bnl.gov/~lebedev/tec/polya.html
      // Gain fluctuation in proportional counters follows Polya distribution.
      // x = G/G_0
      // P(m) = m(m(x)**(m-1)*exp(-m*x)/Gamma(m);
      // original Polya  m = 1.5 (R.Bellazzini and M.A.Spezziga, INFN PI/AE-94/02).
      // Valeri Cherniatin (cherniat@bnlarm.bnl.gov) recomends m=1.38
      // Trs uses x**1.5/exp(x)
      // tss used x**0.5/exp(1.5*x)
      if (St_tpcAltroParamsC::instance()->N(sector - 1) < 0) { // old TPC
        // check if the function has been created
        for (int sec = 1; sec < sector; sec++) {
          if (St_tpcAltroParamsC::instance()->N(sec - 1) < 0 && mShaperResponses[io][sec - 1]) {
            mShaperResponses[io][sector - 1] = mShaperResponses[io][sec - 1];
            break;
          }
        }

        if (! mShaperResponses[io][sector - 1]) {
          mShaperResponses[io][sector - 1] = new TF1F(Form("ShaperFunc_%s_S%02i", Names[io], sector),
              StTpcRSMaker::shapeEI3_I, timeBinMin, timeBinMax, 7);
          mShaperResponses[io][sector - 1]->SetParameters(params3);
          mShaperResponses[io][sector - 1]->SetParNames("t0", "tauF", "tauP", "tauI", "width", "norm", "io");
          mShaperResponses[io][sector - 1]->SetTitle(mShaperResponses[io][sector - 1]->GetName());
          mShaperResponses[io][sector - 1]->GetXaxis()->SetTitle("time (buckets)");
          mShaperResponses[io][sector - 1]->GetYaxis()->SetTitle("signal");
          // Cut tails
          double t = timeBinMax;
          double ymax = mShaperResponses[io][sector - 1]->Eval(0.5);

          for (; t > 5; t -= 1) {
            double r = mShaperResponses[io][sector - 1]->Eval(t) / ymax;

            if (r > 1e-2) break;
          }

          mShaperResponses[io][sector - 1]->SetRange(timeBinMin, t);
          mShaperResponses[io][sector - 1]->Save(timeBinMin, t, 0, 0, 0, 0);
        }

        continue;
      }

      //Altro
      for (int sec = 1; sec < sector; sec++) {
        if (St_tpcAltroParamsC::instance()->N(sec - 1) >= 0 && mShaperResponses[io][sec - 1] ) {
          mShaperResponses[io][sector - 1] = mShaperResponses[io][sec - 1];
          break;
        }
      }

      if (! mShaperResponses[io][sector - 1]) {
        mShaperResponses[io][sector - 1] = new TF1F(Form("ShaperFunc_%s_S%02i", Names[io], sector),
            StTpcRSMaker::shapeEI_I, timeBinMin, timeBinMax, 5);
        mShaperResponses[io][sector - 1]->SetParameters(params0);
        mShaperResponses[io][sector - 1]->SetParNames("t0", "tauI", "width", "norm", "io");
        mShaperResponses[io][sector - 1]->SetTitle(mShaperResponses[io][sector - 1]->GetName());
        mShaperResponses[io][sector - 1]->GetXaxis()->SetTitle("time (buckets)");
        mShaperResponses[io][sector - 1]->GetYaxis()->SetTitle("signal");
        // Cut tails
        double t = timeBinMax;
        double ymax = mShaperResponses[io][sector - 1]->Eval(0.5);

        for (; t > 5; t -= 1) {
          double r = mShaperResponses[io][sector - 1]->Eval(t) / ymax;

          if (r > 1e-2) break;
        }

        mShaperResponses[io][sector - 1]->SetRange(timeBinMin, t);
        mShaperResponses[io][sector - 1]->Save(timeBinMin, t, 0, 0, 0, 0);
      }
    }
  }

  // tss
  mGG = new TF1F("GaitingGridTransperency", "TMath::Max(0.,1-6.27594134307865925e+00*TMath::Exp(-2.87987e-01*(x-1.46222e+01)))", 10, 56);

  if (Debug()) Print();

  memset (hist, 0, sizeof(hist));
  memset (checkList, 0, sizeof(checkList));

  ClusterProfile = true;

  mHeed = fEc(St_TpcResponseSimulatorC::instance()->W());

  if (!ClusterProfile) return;

  int color = 1;
  struct Name_t {
    const char* Name;
    const char* Title;
  };
  const Name_t InOut[6] = {
    {"Inner", "Inner old electronics or iTPC"},
    {"Outer", "Outer old electronics or ITPC"},
    {"InnerX", "Inner new electronics"},
    {"OuterX", "Outer new electronics"},
    {"I", "Inner"},
    {"O", "Outer"}
  };
  const Name_t PadTime[3] = {
    {"Pad", "Pad"},
    {"Time", "Time"},
    {"Row", "Row"},
  };

  for (int io = 2; io < 4; io++) {
    for (int pt = 0; pt < 2; pt++) {
      TString Name(InOut[io].Name); Name += PadTime[pt].Name; Name += "Mc";
      TString Title(InOut[io].Title); Title += PadTime[pt].Title; Title += "Mc";
      hist[io][pt] = (TProfile2D*) gDirectory->Get(Name);

      if (! hist[io][pt]) {
        hist[io][pt] = new TProfile2D(Name, Title, nx[pt], xmin[pt], xmax[pt], nz, zmin, zmax, "");
        hist[io][pt]->SetMarkerStyle(20);
        hist[io][pt]->SetMarkerColor(color++);
      }
    }
  }

  hist[4][0] = new TProfile2D("dEdxCorSecRow", "dEdx correction versus sector and row",
                              max_sectors, 0.5, max_sectors + 0.5,
                              St_tpcPadConfigC::instance()->numberOfRows(20), 0.5, St_tpcPadConfigC::instance()->numberOfRows(20) + 0.5, "");
  hist[4][1] = new TProfile2D("GainSecRow", "Overall gain versus sector and row",
                              max_sectors, 0.5, max_sectors + 0.5,
                              St_tpcPadConfigC::instance()->numberOfRows(20), 0.5, St_tpcPadConfigC::instance()->numberOfRows(20) + 0.5, "");
  const Name_t Checks[21] = {
    {"dEGeant", "dE in Geant"}, // 0
    {"dSGeant", "ds in Geant"}, // 1
    {"Gain", "Gas Gain after Voltage"}, // 2
    {"GainMc", "Gas Gain after MC correction"}, // 3
    {"dEdxCor", "correction of dEdx"}, // 4
    {"lgam", "lgam"}, // 5
    {"NPGEANT", "no. of primary electros from GEANT"}, // 6
    {"NP", "no. of primary electros"}, // 7
    {"Nt", "total no. of electors per cluster"}, // 8
    {"Qav", "Gas gain flactuations"}, // 9
    {"localYDirectionCoupling", "localYDirectionCoupling"}, //10
    {"n0", "No. electrons per primary interaction"}, //11
    {"padGain", "padGain"}, // 12
    {"localXDirectionCoupling", "localXDirectionCoupling"}, // 13
    {"XYcoupling", "XYcoupling"}, //14
    {"dE", "dE"}, // 15
    {"dS", "dS"}, // 16
    {"adc", "adc"}, // 17
    {"NE", "Total no. of generated electors"}, // 18
    {"dECl", "Total log(signal/Nt) in a cluster versus Wire Index"}, // 19
    {"nPdT", "log(Total no. of conducting electrons) - log(no. of primary one) versus log(no. primary electrons)"} // 20
  };
  const int Npbins  = 151;
  const int NpbinsL =  10;
  const double Xmax = 1e5;
  double    dX = TMath::Log(Xmax / 10) / (Npbins - NpbinsL);
  double* pbins = new double[Npbins];
  double* pbinsL =  new double[Npbins];
  pbins[0] = 0.5;
  pbinsL[0] = TMath::Log(pbins[0]);

  for (int bin = 1; bin < Npbins; bin++) {
    if (bin <= NpbinsL) {
      pbins[bin] = pbins[bin - 1] + 1;
    }
    else if (bin == Npbins - 1) {
      pbins[bin] = 1e5;
    }
    else {
      int nM = 0.5 * (pbins[NpbinsL - 2] + pbins[NpbinsL - 1]) * TMath::Exp(dX * (bin - NpbinsL));
      double dbin = tpcrs::irint(nM - pbins[bin - 1]);

      if (dbin < 1.0) dbin = 1.0;

      pbins[bin] = pbins[bin - 1] + dbin;
    }

    pbinsL[bin] = TMath::Log(pbins[bin]);
  }

  for (int io = 0; io < 2; io++) {
    for (int i = 0; i < nChecks; i++) {
      TString Name(Checks[i].Name); Name += InOut[4 + io].Name;
      TString Title(Checks[i].Title); Title += InOut[4 + io].Title;

      if      (i == 11) checkList[io][i] = new TH2D(Name, Title, nz, zmin, zmax, 100, -0.5, 99.5);
      else if (i == 19) checkList[io][i] = new TH2D(Name, Title, 173, -.5, 172.5, 200, -10, 10);
      //      else if (i == 20) checkList[io][i] = new TH2D(Name,Title,Npbins-1,pbinsL,Npbins-1,pbinsL);
      else if (i == 20) checkList[io][i] = new TH2D(Name, Title, Npbins - 1, pbinsL, 500, -2.0, 8.0);
      else              checkList[io][i] = new TProfile(Name, Title, nz, zmin, zmax, "");
    }
  }

  delete [] pbins;
  delete [] pbinsL;
}


inline bool operator< (const g2t_tpc_hit_st& lhs, const g2t_tpc_hit_st& rhs)
{
  // sectors
  if ((lhs.volume_id % 100000) / 100 != (rhs.volume_id % 100000) / 100)
    return (lhs.volume_id % 100000) / 100 < (rhs.volume_id % 100000) / 100;

  // track id
  if (lhs.track_p != rhs.track_p)
    return lhs.track_p < rhs.track_p;

  // pad rows
  //  if (lhs.volume_id%100 != rhs.volume_id%100) return lhs.volume_id%100 - rhs.volume_id%100;
  // track length
  return lhs.length < rhs.length;
}


void StTpcRSMaker::Make(const std::vector<g2t_tpc_hit_st>& g2t_tpc_hit,
                        std::vector<g2t_track_st>& g2t_track,
                        const std::vector<g2t_vertex_st>& g2t_vertex, tpcrs::DigiData& digi_data)
{
  static int nCalls = 0;
  gRandom->SetSeed(2345 + nCalls++);

  static int max_sectors = 24;
  // constants
  static int iBreak = 0;
#ifdef __DEBUG__
  if (Debug() % 10) {
    gBenchmark->Reset();
    gBenchmark->Start("TpcRS");
    LOG_INFO << "\n -- Begin TpcRS Processing -- \n";
  }
#endif
  double vminI = St_tpcGainCorrectionC::instance()->Struct(1)->min;
  double vminO = St_tpcGainCorrectionC::instance()->Struct(0)->min;

  // TODO: Confirm proper handling of empty input containers
  int no_tpc_hits = g2t_tpc_hit.size();

  int NoTpcTracks = g2t_track.size();

  std::vector<int> mNoTpcHitsAll(NoTpcTracks + 1);
  std::vector<int> mNoTpcHitsReal(NoTpcTracks + 1);

  g2t_track_st* geantTrack = &g2t_track.front();

  const g2t_vertex_st* geantVertex = &g2t_vertex.front();
  int NV = g2t_vertex.size();

  g2t_tpc_hit_st* tpc_hit_begin = const_cast<g2t_tpc_hit_st*>(&g2t_tpc_hit.front());
  g2t_tpc_hit_st* tpc_hit = tpc_hit_begin;

  if (m_TpcdEdxCorrection) {
    St_tpcGainCorrectionC::instance()->Struct(0)->min = -500;
    St_tpcGainCorrectionC::instance()->Struct(1)->min = -500;

    if (Debug()) {
      LOG_INFO << "Reset min for gain Correction to I/O\t"
               << St_tpcGainCorrectionC::instance()->Struct(1)->min
               << "\t"
               << St_tpcGainCorrectionC::instance()->Struct(0)->min
               << " (V)\n";
    }
  }

  mNSplittedHits = 0;
  // sort
  // initialize original index locations
  std::vector<size_t> sorted_index(no_tpc_hits);
  std::iota(sorted_index.begin(), sorted_index.end(), 0);
  std::stable_sort(sorted_index.begin(), sorted_index.end(), [&g2t_tpc_hit](size_t i1, size_t i2) {return g2t_tpc_hit[i1] < g2t_tpc_hit[i2];}); 

  int sortedIndex = 0;

  for (int sector = 1; sector <= max_sectors; sector++) {
    int nHitsInTheSector = 0;
    free(m_SignalSum); m_SignalSum = 0;
    ResetSignalSum(sector);

    // it is assumed that hit are ordered by sector, trackId, pad rows, and track length
    for (; sortedIndex < no_tpc_hits; sortedIndex++) {
      int indx = sorted_index[sortedIndex];

      if (indx < 0) break;

      tpc_hit = tpc_hit_begin + indx;
      int volId = tpc_hit->volume_id % 10000;
      int iSector = volId / 100;

      if (iSector != sector) {
        if (! ( iSector > sector ) ) {
          LOG_ERROR << "StTpcRSMaker::Make: g2t_tpc_hit table has not been ordered by sector no. " << sector << '\n';
          assert( iSector > sector );
        }

        break;
      }

      if (tpc_hit->volume_id <= 0 || tpc_hit->volume_id > 1000000) continue;

      int Id  = tpc_hit->track_p;
      int id3 = 0, ipart = 8, charge = 1;
      double mass = 0;

      if (geantTrack) {
        id3        = geantTrack[Id - 1].start_vertex_p;
        assert(id3 > 0 && id3 <= NV);
        ipart      = geantTrack[Id - 1].ge_pid;
        charge     = (int) geantTrack[Id - 1].charge;
        StParticleDefinition* particle = StParticleTable::instance()->findParticleByGeantId(ipart);

        if (particle) {
          mass = particle->mass();
          charge = particle->charge();
        }
      }

      if (ipart == Laserino || ipart == Chasrino) {
        charge = 0;
      }
      else {
        if (ipart == 1) {// gamma => electron
          ipart = 3;
          charge = -1;
        }

        if (charge == 0) {
          continue;
        }
      } // special treatment for electron/positron

      if (ipart == 2) charge =  101;

      if (ipart == 3) charge = -101;

      // Track segment to propagate
      enum {NoMaxTrackSegmentHits = 100};
      static HitPoint_t TrackSegmentHits[NoMaxTrackSegmentHits];
      msMin = 9999;
      msMax = -9999;
      int nSegHits = 0;
      int sIndex = sortedIndex;

      if (Debug() > 13) LOG_INFO << "sortedIndex = " << sortedIndex << "\tno_tpc_hits = " << no_tpc_hits << '\n';

      int ID = 0;
      int TrackDirection = 0; // 0 - increase no of row, 1 - decrease no of. row.

      for (nSegHits = 0, sIndex = sortedIndex; sIndex < no_tpc_hits && nSegHits < NoMaxTrackSegmentHits; sIndex++)
      {
        indx = sorted_index[sIndex];
        g2t_tpc_hit_st* tpc_hitC = tpc_hit_begin + indx;

        if ((tpc_hitC->volume_id % 10000) / 100 != sector) break;

        if (ID > 0 && ID != tpc_hitC->track_p) break;

        ID = tpc_hitC->track_p;

        if (nSegHits == 1) { // No Loopers !
          if (TrackSegmentHits[nSegHits - 1].tpc_hitC->volume_id % 100 <= tpc_hitC->volume_id % 100) {
            TrackDirection = 0;
          }
          else {
            TrackDirection = 1;
          }
        }
        else if (nSegHits > 1) {
          if ((! TrackDirection && TrackSegmentHits[nSegHits - 1].tpc_hitC->volume_id % 100 > tpc_hitC->volume_id % 100) ||
              (  TrackDirection && TrackSegmentHits[nSegHits - 1].tpc_hitC->volume_id % 100 < tpc_hitC->volume_id % 100))
            break;
        }

        if (Debug() > 13) LOG_INFO << "sIndex = " << sIndex << "\tindx = " << indx << "\ttpc_hitC = " << tpc_hitC << '\n';

        TrackSegmentHits[nSegHits].indx = indx;
        TrackSegmentHits[nSegHits].s = tpc_hitC->length;

        if (tpc_hitC->length == 0 && nSegHits > 0) {
          TrackSegmentHits[nSegHits].s = TrackSegmentHits[nSegHits - 1].s + TrackSegmentHits[nSegHits].tpc_hitC->ds;
        }

        // Increase hit counters
        mNoTpcHitsAll[tpc_hitC->track_p - 1]++;
        // Account hits which can be splitted
        if (tpc_hitC->volume_id < 10000) mNoTpcHitsReal[tpc_hitC->track_p - 1]++;

        TrackSegment2Propagate(tpc_hitC, &geantVertex[id3 - 1], TrackSegmentHits[nSegHits]);

        if (TrackSegmentHits[nSegHits].Pad.timeBucket() < 0 || TrackSegmentHits[nSegHits].Pad.timeBucket() > max_timebins) continue;

        nSegHits++;
      }

      if (! nSegHits) continue;

      if (Debug() >= 10) {
        PrPP(Make, nSegHits);

        for (int s = 0; s < nSegHits; s++) {
          LOG_INFO << "Seg[" << Form("%2i", s) << "]\tId " << TrackSegmentHits[s].TrackId << "\ts = " << TrackSegmentHits[s].s
               << "\tvolumeID :" <<  Form("%6i", TrackSegmentHits[s].tpc_hitC->volume_id) << "\t" << TrackSegmentHits[s].Pad
               << "\ts1/s2 = " << TrackSegmentHits[s].tpc_hitC->length - TrackSegmentHits[s].tpc_hitC->ds / 2
               << "\t" << TrackSegmentHits[s].tpc_hitC->length + TrackSegmentHits[s].tpc_hitC->ds / 2 << "\tds = " << TrackSegmentHits[s].tpc_hitC->ds
               << '\n';
        }
      }

      sortedIndex = sIndex - 1; // Irakli 05/06/19, reduce extra step in for loop
      double s = msMin;
      memset (rowsdE, 0, sizeof(rowsdE));

      for (int iSegHits = 0; iSegHits < nSegHits && s < msMax; iSegHits++) {
        memset (rowsdEH, 0, sizeof(rowsdEH));
        g2t_tpc_hit_st* tpc_hitC = TrackSegmentHits[iSegHits].tpc_hitC;
        tpc_hitC->adc = 0;
        volId = tpc_hitC->volume_id % 100000;
        int row = TrackSegmentHits[iSegHits].coorLS.fromRow();
        int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;
        // switch between Inner / Outer Sector paramters
        // Extra correction for simulation with respect to data
        int iowe = 0;

        if (sector  > 12) iowe += 4;

        if (io) iowe += 2;

        float*  AdditionalMcCorrection = St_TpcResponseSimulatorC::instance()->SecRowCor();
        float*  AddSigmaMcCorrection   = St_TpcResponseSimulatorC::instance()->SecRowSig();
        // Generate signal
        double sigmaJitterT     = St_TpcResponseSimulatorC::instance()->SigmaJitterTI();
        double sigmaJitterX     = St_TpcResponseSimulatorC::instance()->SigmaJitterXI();

        if (io) { // Outer
          sigmaJitterT            = St_TpcResponseSimulatorC::instance()->SigmaJitterTO();
          sigmaJitterX            = St_TpcResponseSimulatorC::instance()->SigmaJitterXO();
        }

        // Generate signal
        double Gain = St_tss_tssparC::instance()->gain(sector, row);
        mShaperResponse = mShaperResponses[io][sector - 1];

        if (ClusterProfile) {
          checkList[io][2]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), Gain);
        }

        double GainXCorrectionL = AdditionalMcCorrection[iowe] + row * AdditionalMcCorrection[iowe + 1];
        Gain *= TMath::Exp(-GainXCorrectionL);
        double GainXSigma = AddSigmaMcCorrection[iowe] + row * AddSigmaMcCorrection[iowe + 1];

        if (GainXSigma > 0) Gain *= TMath::Exp(gRandom->Gaus(0., GainXSigma));

        if (ClusterProfile) {
          checkList[io][3]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), Gain);
        }

        // dE/dx correction
        double dEdxCor = dEdxCorrection(TrackSegmentHits[iSegHits]);
#ifdef __DEBUG__
        if (TMath::IsNaN(dEdxCor)) {
          iBreak++;
        }
#endif
        if (dEdxCor <= 0.) continue;

        if (ClusterProfile) {
          checkList[io][4]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), dEdxCor);
          hist[4][0]->Fill(TrackSegmentHits[iSegHits].Pad.sector(), TrackSegmentHits[iSegHits].Pad.row(), dEdxCor);
        }

        // Apply Gating Grid
        if (TrackSegmentHits[iSegHits].Pad.timeBucket() > mGG->GetXmin() &&
            TrackSegmentHits[iSegHits].Pad.timeBucket() < mGG->GetXmax()) {
          dEdxCor *= mGG->Eval(TrackSegmentHits[iSegHits].Pad.timeBucket());
        }

        if (dEdxCor < minSignal) continue;

        // Initialize propagation
        float BField[3] = {(float ) TrackSegmentHits[iSegHits].BLS.position().x(),
                             (float ) TrackSegmentHits[iSegHits].BLS.position().y(),
                             (float ) TrackSegmentHits[iSegHits].BLS.position().z()
                            };
        // Magnetic field BField must be in kilogauss
        // kilogauss = 1e-1*tesla = 1e-1*(volt*second/meter2) = 1e-1*(1e-6*1e-3*1/1e4) = 1e-14
        StPhysicalHelix track(TrackSegmentHits[iSegHits].dirLS.position(),
                              TrackSegmentHits[iSegHits].coorLS.position(),
                              BField[2]* 1e-14 * charge, 1);
        StThreeVector<double> unit = TrackSegmentHits[iSegHits].dirLS.position().unit();
        double* cxyz = unit.xyz();
        double L2L[9] = {
                     cxyz[2], - cxyz[0]*cxyz[2], cxyz[0],
                     cxyz[0], - cxyz[1]*cxyz[2], cxyz[1],
                     0.0,   cxyz[0]*cxyz[0] + cxyz[1]*cxyz[1], cxyz[2]};
#ifdef __DEBUG__
        if (Debug() > 11) PrPP(Make, track);
#endif
        double dStep =  TMath::Abs(tpc_hitC->ds);
        double s_low   = -dStep / 2;
        double s_upper = s_low + dStep;
        double newPosition = s_low;
        static StThreeVector<double> normal(0, 1, 0);
        static StTpcCoordinateTransform transform(gStTpcDb);
        StThreeVector<double> rowPlane(0, transform.yFromRow(TrackSegmentHits[iSegHits].Pad.sector(), TrackSegmentHits[iSegHits].Pad.row()), 0);
        double sR = track.pathLength(rowPlane, normal);

        if (sR < 1e10) {
          PrPP(Maker, sR);
          PrPP(Make, TrackSegmentHits[iSegHits].coorLS);
          StThreeVector<double> xyzP = track.at(sR);
          TrackSegmentHits[iSegHits].coorLS.setPosition(xyzP);
          PrPP(Make, TrackSegmentHits[iSegHits].coorLS);
          // don't useT0, don't useTau
          PrPP(Make, TrackSegmentHits[iSegHits].Pad);
          transform(TrackSegmentHits[iSegHits].coorLS, TrackSegmentHits[iSegHits].Pad, false, false); // don't use T0, don't use Tau
          PrPP(Make, TrackSegmentHits[iSegHits].Pad);
        }

        int ioH = io;

        if (St_tpcAltroParamsC::instance()->N(sector - 1) >= 0) ioH += 2;

        TotalSignal  = 0;
        double lgam = tpc_hitC->lgam;

        if (ClusterProfile) {
          checkList[io][5]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), lgam);
        }

        double gamma = TMath::Power(10., lgam) + 1;
        double betaGamma = TMath::Sqrt(gamma * gamma - 1.);
        StThreeVector<double>       pxyzG(tpc_hitC->p[0], tpc_hitC->p[1], tpc_hitC->p[2]);
        double bg = 0;
        static const double m_e = .51099907e-3;
        static const double eV = 1e-9; // electronvolt in GeV
        double eKin = -1;
#ifdef __STOPPED_ELECTRONS__

        if (mass > 0) {
          bg = pxyzG.mag() / mass;

          // special case stopped electrons
          if (tpc_hitC->ds < 0.0050 && tpc_hitC->de < 0) {
            int Id         = tpc_hitC->track_p;
            int ipart      = geantTrack[Id - 1].ge_pid;

            if (ipart == 3) {
              eKin = -tpc_hitC->de;
              gamma = eKin / m_e + 1;
              bg = TMath::Sqrt(gamma * gamma - 1.);
            }
          }
        }

#else /* ! __STOPPED_ELECTRONS__ */

        if (mass > 0) bg = pxyzG.mag() / mass;

#endif /* __STOPPED_ELECTRONS__ */

        if (bg > betaGamma) betaGamma = bg;

        double bg2 = betaGamma * betaGamma;
        gamma = TMath::Sqrt(bg2 + 1.);
        double Tmax;

        if (mass < 2 * m_e) {
          if (charge > 0) Tmax =     m_e * (gamma - 1);
          else            Tmax = 0.5 * m_e * (gamma - 1);
        }
        else {
          double r = m_e / mass;
          Tmax = 2 * m_e * bg2 / (1 + 2 * gamma * r + r * r);
        }

        if (Tmax > mCutEle) Tmax = mCutEle;

        double padH = TrackSegmentHits[iSegHits].Pad.pad();
        double tbkH = TrackSegmentHits[iSegHits].Pad.timeBucket();
        tpc_hitC->pad = padH;
        tpc_hitC->timebucket = tbkH;
        pad0 = tpcrs::irint(padH + xmin[0]);
        tbk0 = tpcrs::irint(tbkH + xmin[1]);
        double OmegaTau = St_TpcResponseSimulatorC::instance()->OmegaTau() *
                            TrackSegmentHits[iSegHits].BLS.position().z() / 5.0; // from diffusion 586 um / 106 um at B = 0/ 5kG
        double NP = TMath::Abs(tpc_hitC->de) / (St_TpcResponseSimulatorC::instance()->W() * eV *
                      St_TpcResponseSimulatorC::instance()->Cluster()); // from GEANT

        if (ClusterProfile) {
          checkList[io][6]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), NP);
        }

        double driftLength = TMath::Abs(TrackSegmentHits[iSegHits].coorLS.position().z());
        double D = 1. + OmegaTau * OmegaTau;
        double SigmaT = St_TpcResponseSimulatorC::instance()->transverseDiffusion() *  TMath::Sqrt(   driftLength / D);

        //	double SigmaL = St_TpcResponseSimulatorC::instance()->longitudinalDiffusion()*TMath::Sqrt(2*driftLength  );
        if (sigmaJitterX > 0) {SigmaT = TMath::Sqrt(SigmaT * SigmaT + sigmaJitterX * sigmaJitterX);}

        double SigmaL = St_TpcResponseSimulatorC::instance()->longitudinalDiffusion() * TMath::Sqrt(   driftLength  );
        double NoElPerAdc = St_TpcResponseSimulatorC::instance()->NoElPerAdc();

        if (NoElPerAdc <= 0) {
          if (St_tpcPadConfigC::instance()->iTPC(sector) && St_tpcPadConfigC::instance()->IsRowInner(sector, row)) {
            NoElPerAdc = St_TpcResponseSimulatorC::instance()->NoElPerAdcX(); // iTPC
          }
          else if (St_tpcPadConfigC::instance()->IsRowInner(sector, row)) {
            NoElPerAdc = St_TpcResponseSimulatorC::instance()->NoElPerAdcI(); // inner TPX
          }
          else {
            NoElPerAdc = St_TpcResponseSimulatorC::instance()->NoElPerAdcO(); // outer TPX
          }
        }

#ifndef __NO_1STROWCORRECTION__
        if (row == 1) dEdxCor *= TMath::Exp(St_TpcResponseSimulatorC::instance()->FirstRowC());
#endif
        mGainLocal = Gain / dEdxCor / NoElPerAdc; // Account dE/dx calibration
        // end of dE/dx correction
        // generate electrons: No. of primary clusters per cm

        if (mdNdx || mdNdxL10) {
          NP = GetNoPrimaryClusters(betaGamma, charge); // per cm
#ifdef __DEBUG__
          if (NP <= 0.0) {
            iBreak++; continue;
          }
#endif
        }

        if (ClusterProfile) {
          checkList[io][7]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), NP);
        }

        int nP = 0;
        double dESum = 0;
        double dSSum = 0;
        int   nTotal = 0;
        memset (padsdE, 0, sizeof(padsdE));
        memset (tbksdE,  0, sizeof(tbksdE));
        float dEr = 0;
        std::vector<float> rs(10);

        do {// Clusters
          float dS = 0;
          float dE = 0;
          static double cLog10 = TMath::Log(10.);

          if (eKin >= 0.0) {
            if (eKin == 0.0) break;

            gamma = eKin / m_e + 1;
            bg = TMath::Sqrt(gamma * gamma - 1.);
            Tmax = 0.5 * m_e * (gamma - 1);

            if (Tmax <= St_TpcResponseSimulatorC::instance()->W() / 2 * eV) break;

            NP = GetNoPrimaryClusters(betaGamma, charge);
            dE = TMath::Exp(cLog10 * mdNdEL10->GetRandom());
          }
          else {
            if (charge) {
              dS = - TMath::Log(gRandom->Rndm()) / NP;

              if (mdNdEL10) dE = TMath::Exp(cLog10 * mdNdEL10->GetRandom());
              else          dE = St_TpcResponseSimulatorC::instance()->W() *
                                   gRandom->Poisson(St_TpcResponseSimulatorC::instance()->Cluster());
            }
            else { // charge == 0 geantino
              // for Laserino assume dE/dx = 25 keV/cm;
              dE = 10; // eV
              dS = dE * eV / (TMath::Abs(mLaserScale * tpc_hitC->de / tpc_hitC->ds));
            }
          }

#ifdef __DEBUG__
          if (Debug() > 12) {
            LOG_INFO << "s_low/s_upper/dSD\t" << s_low << "/\t" << s_upper << "\t" << dS <<  '\n';
          }
#endif
          double E = dE * eV;
          newPosition += dS;

          if (newPosition > s_upper) break;

          if (dE < St_TpcResponseSimulatorC::instance()->W() / 2 || E > Tmax) continue;

          if (eKin > 0) {
            if (eKin >= E) {eKin -= E;}
            else {E = eKin; eKin = 0; dE = E / eV;}
          }

          dESum += dE;
          dSSum += dS;
          nP++;
#ifdef __DEBUG__
          if (Debug() > 12) {
            LOG_INFO << "dESum = " << dESum << " /\tdSSum " << dSSum << " /\t newPosition " << newPosition << '\n';
          }
#endif
          double xRange = 0;

          if (dE > ElectronRangeEnergy) xRange = ElectronRange * TMath::Power((dE + dEr) / ElectronRangeEnergy, ElectronRangePower);

          int Nt = 0; // HeedCsize(dE, dEr,rs);
          float dET = dE + dEr;
          dEr = dET;
          float EC;
          int   Nr = 0;

          if (xRange > 0) {Nr = rs.size();}

          while ((EC = mHeed->GetRandom()) < dEr) {
            dEr -= EC;

            if (Nr) {
              if (Nr <= Nt) { Nr = 2 * Nt + 1; rs.resize(Nr); }
              rs[Nt] = 1 - dEr / dET;
            }

            Nt++;
          }

          if (!Nt) continue;

          if (ClusterProfile) {
            checkList[io][8]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), Nt);
            checkList[io][11]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), Nt);
          }

          StThreeVector<double> xyzC = track.at(newPosition);
          double phiXY = 2 * TMath::Pi() * gRandom->Rndm();
          double rX = TMath::Cos(phiXY);
          double rY = TMath::Sin(phiXY);
          double sigmaT = SigmaT;
          double sigmaL = SigmaL;
          TotalSignalInCluster = 0;
          int WireIndex = 0;

          for (int ie = 0; ie < Nt; ie++) {
            nTotal++;
            QAv = mPolya[io]->GetRandom();
            // transport to wire
            gRandom->Rannor(rX, rY);
            StTpcLocalSectorCoordinate xyzE(xyzC.x() + rX * sigmaT,
                                            xyzC.y() + rY * sigmaT,
                                            xyzC.z() + gRandom->Gaus(0, sigmaL), sector, row);

            if (xRange > 0) {
              double xR = rs[ie] * xRange;
              double xyzRangeL[3] = {xR * rX, xR * rY, 0.};
              double xyzR[3] = {0};
              TCL::mxmpy(L2L, xyzRangeL, xyzR, 3, 3, 1);
#ifdef __DEBUG__
              if (Debug() > 13) {
                LOG_INFO << "xyzRangeL: " << xyzRangeL[0] << " " << xyzRangeL[1] << " " << xyzRangeL[2] << '\n';
                LOG_INFO << "L2L: " << L2L[0] << " " << L2L[1] << " " << L2L[2]
                                    << L2L[3] << " " << L2L[4] << " " << L2L[3]
                                    << L2L[6] << " " << L2L[7] << " " << L2L[8] << '\n';
                LOG_INFO << "xyzR: " << xyzR[0] << " " << xyzR[1] << " " << xyzR[2] << '\n';
              }
#endif
              for (int i=0; i<3; i++) xyzE.position().xyz()[i] += xyzR[i];
            }

            double y = xyzE.position().y();
            double alphaVariation = InnerAlphaVariation[sector - 1];

            // Transport to wire
            if (y <= lastInnerSectorAnodeWire) {
              WireIndex = tpcrs::irint((y - firstInnerSectorAnodeWire) / anodeWirePitch) + 1;
#ifndef __NO_1STROWCORRECTION__
              if (St_tpcPadConfigC::instance()->iTPC(sector)) {// two first and two last wires are removed, and 3rd wire is fat wiere
                if (WireIndex <= 3 || WireIndex >= numberOfInnerSectorAnodeWires - 3) continue;
              }
              else {   // old TPC the first and last wires are fat ones
                if (WireIndex <= 1 || WireIndex >= numberOfInnerSectorAnodeWires) continue;
              }
#else
              if (WireIndex <= 1 || WireIndex >= numberOfInnerSectorAnodeWires) continue; // to check the 1-st pad row effect
#endif
              yOnWire = firstInnerSectorAnodeWire + (WireIndex - 1) * anodeWirePitch;
            }
            else {   // the first and last wires are fat ones
              WireIndex = tpcrs::irint((y - firstOuterSectorAnodeWire) / anodeWirePitch) + 1;

              if (WireIndex <= 1 || WireIndex >= numberOfOuterSectorAnodeWires) continue;

              yOnWire = firstOuterSectorAnodeWire + (WireIndex - 1) * anodeWirePitch;
              alphaVariation = OuterAlphaVariation[sector - 1];
            }

            double distanceToWire = y - yOnWire; // Calculated effective distance to wire affected by Lorentz shift
            xOnWire = xyzE.position().x();
            zOnWire = xyzE.position().z();
            // Grid plane (1 mm spacing) focusing effect + Lorentz angle in drift volume
            int iGridWire = (int ) TMath::Abs(10.*distanceToWire);
            double dist2Grid = TMath::Sign(0.05 + 0.1 * iGridWire, distanceToWire); // [cm]
            // Ground plane (1 mm spacing) focusing effect
            int iGroundWire = (int ) TMath::Abs(10.*dist2Grid);
            double distFocused = TMath::Sign(0.05 + 0.1 * iGroundWire, dist2Grid);
            // OmegaTau near wires taken from comparison with data
            double tanLorentz = OmegaTau / St_TpcResponseSimulatorC::instance()->OmegaTauScaleO();

            if (y < firstOuterSectorAnodeWire) tanLorentz = OmegaTau / St_TpcResponseSimulatorC::instance()->OmegaTauScaleI();

            xOnWire += distFocused * tanLorentz; // tanLorentz near wires taken from comparison with data
            zOnWire += TMath::Abs(distFocused);

            if (! iGroundWire ) QAv *= TMath::Exp( alphaVariation);
            else                QAv *= TMath::Exp(-alphaVariation);

            if (ClusterProfile) {
              checkList[io][9]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), QAv);
            }

            double dY    = mChargeFraction[io][sector - 1]->GetXmax();
            double yLmin = yOnWire - dY;
            double yLmax = yLmin + 2 * dY;
            int    rowMin  = transform.rowFromLocalY(yLmin, sector);
            int    rowMax  = transform.rowFromLocalY(yLmax, sector);
            double yRmin = transform.yFromRow(sector, rowMin) - St_tpcPadConfigC::instance()->PadLengthAtRow(sector, rowMin) / 2;
            double yRmax = transform.yFromRow(sector, rowMax) + St_tpcPadConfigC::instance()->PadLengthAtRow(sector, rowMax) / 2;

            if (yRmin > yLmax || yRmax < yLmin) {
              iBreak++; continue;
            }

            GenerateSignal(TrackSegmentHits[iSegHits], sector, rowMin, rowMax, sigmaJitterT, sigmaJitterX);
          }  // electrons in Cluster

          if (ClusterProfile) {
            if (TotalSignalInCluster > 0 && checkList[io][19]) {
              checkList[io][19]->Fill(WireIndex, TMath::Log(TotalSignalInCluster));

            }
          }

          TotalSignal += TotalSignalInCluster;
        }
        while (true);   // Clusters

        //	tpc_hitC->adc = -99;
        if (dESum > 0 && dSSum) {
#ifdef __DEBUG__
          if (Debug() > 12) {
            LOG_INFO << "sIndex = " << sIndex << " volId = " << volId
                     << " dESum = " << dESum << " /\tdSSum " << dSSum << " /\t TotalSignal " << TotalSignal << '\n';
          }
#endif
          tpc_hitC->de = dESum * eV;
          tpc_hitC->ds = dSSum;
          //	  tpc_hitC->adc = TotalSignal;
          tpc_hitC->np = nP;

          if (ClusterProfile) {
            if (TotalSignal > 0) {
              if (hist[ioH][0]) {
                for (int p = 0; p < kPadMax; p++)
                  hist[ioH][0]->Fill((p + pad0) - padH, TrackSegmentHits[iSegHits].xyzG.position().z(), padsdE[p] / TotalSignal);
              }

              if (hist[ioH][1]) {
                for (int t = 0; t < kTimeBacketMax; t++)
                  hist[ioH][1]->Fill((t + tbk0 + 0.5) - tbkH, TrackSegmentHits[iSegHits].xyzG.position().z(), tbksdE[t] / TotalSignal);
              }
            }

            checkList[io][15]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), tpc_hitC->de);
            checkList[io][16]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), tpc_hitC->ds);
            checkList[io][18]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), nTotal);

            if (nP > 0 && nTotal > 0)
              checkList[io][20]->Fill(TMath::Log(nP), TMath::Log(nTotal) - TMath::Log(nP));
          }
        }

        nHitsInTheSector++;
      } // end do loop over segments for a given particle

      for (int iSegHits = 0; iSegHits < nSegHits; iSegHits++) {
        g2t_tpc_hit_st* tpc_hitC = TrackSegmentHits[iSegHits].tpc_hitC;

        if (tpc_hitC->volume_id > 10000) continue;

        int row = tpc_hitC->volume_id % 100;
        tpc_hitC->adc += rowsdE[row - 1];
        int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;

        if (checkList[io][17])
          checkList[io][17]->Fill(TrackSegmentHits[iSegHits].xyzG.position().z(), tpc_hitC->adc);
      }
    }  // hits in the sector

    if (nHitsInTheSector) {
      DigitizeSector(sector, digi_data);

      if (Debug()) LOG_INFO << "StTpcRSMaker: Done with sector\t" << sector << " total no. of hit = " << nHitsInTheSector << '\n';
    }
  } // sector

  if (m_SignalSum) {free(m_SignalSum); m_SignalSum = 0;}

  if (Debug() % 10) gBenchmark->Show("TpcRS");

  if (m_TpcdEdxCorrection) {
    St_tpcGainCorrectionC::instance()->Struct(1)->min = vminI;
    St_tpcGainCorrectionC::instance()->Struct(0)->min = vminO;

    if (Debug()) {
      LOG_INFO << "Reset min for gain Correction to I/O\t"
               << St_tpcGainCorrectionC::instance()->Struct(1)->min
               << "\t"
               << St_tpcGainCorrectionC::instance()->Struct(0)->min
               << " (V)\n";
    }
  }

    // Reset no. Tpc hits in g2t_track
    geantTrack = &g2t_track.front();

    for (int i = 0; i < NoTpcTracks; i++, geantTrack++) {
      int Id = geantTrack->id;
      geantTrack->n_tpc_hit = (mNoTpcHitsReal[Id - 1] << 8) + (0xff & mNoTpcHitsAll[Id - 1]);
    }
}


double StTpcRSMaker::GetNoPrimaryClusters(double betaGamma, int charge)
{
  double beta = betaGamma / TMath::Sqrt(1.0 + betaGamma * betaGamma);
  double dNdx = 0;

  if      (mdNdx   ) dNdx = mdNdx->Interpolate(betaGamma);
  else if (mdNdxL10) dNdx = mdNdxL10->Interpolate(TMath::Log10(betaGamma));

  double Q_eff = TMath::Abs(charge % 100);

  if (Q_eff > 1)   {
    // Effective charge from GEANT ghion.F
    double w1 = 1.034 - 0.1777 * TMath::Exp(-0.08114 * Q_eff);
    double w2 = beta * TMath::Power(Q_eff, -2. / 3.);
    double w3 = 121.4139 * w2 + 0.0378 * TMath::Sin(190.7165 * w2);
    Q_eff      *= 1. - w1 * TMath::Exp(-w3);
  }

  return Q_eff * Q_eff * dNdx;
}


double StTpcRSMaker::ShaperFunc(double* x, double* par)
{
  double tau = par[0];
  double width = par[1];
  double p = par[2];
  double t = x[0] * width / tau;
  double Delta = width / tau;
  double t1 = t - Delta / 2.;
  double t2 = t1 + Delta;
  double val = TMath::Gamma(p, t2) - TMath::Gamma(p, t1);
  return val;
}


double StTpcRSMaker::PadResponseFunc(double* x, double* par)
{
  double CrossTalk = 0;
  double Value = 0;
  double X = par[5] * x[0];

  if (CrossTalk > 0) {
    for (int i = -1; i <= 1; i++) {
      double xx = X + par[5] * i;

      if (i == 0) Value += (1. - 2.*CrossTalk) * Gatti(&xx, par);
      else        Value +=          CrossTalk * Gatti(&xx, par);
    }
  }
  else   Value = Gatti(&X, par);

  return Value;
}


double StTpcRSMaker::Gatti(double* x, double* par)
{
  /************************************************************************
   *  Function    : generates the cathode signal using                    *
   *                the single-parameter Gatti formula:                   *
   *                              1 - tanh(K2 * lambda)**2                *
   *     GFunc(lambda) = K1 * -------------------------------             *
   *                           1 + K3 * tanh (K2 *lambda)**2              *
   *     lambda = x/h, h is anode cathode spacing                         *
   *                                                                      *
   *     K2 = pi/2*(1 - 0.5*sqrt(K3))                                     *
   *                                                                      *
   *              K2*sqrt(K3)                                             *
   *     K1 = -------------------                                         *
   *            4 * atan(sqrt(K3))                                        *
   *                                                                      *
   *  References  : E.Gatti, A.Longoni, NIM 163 (1979) 82-93.             *
   *  Authors : V.Balagura,V.Cherniatin,A.Chikanian                       *
   ************************************************************************/
  double y = x[0];   // distance to center of strip [cm]
  double w = par[0]; // w = width of pad
  double h = par[1]; // h = Anode-Cathode gap
  double K3  = par[3];
  double lambda = y / h;
  double K2 = TMath::PiOver2() * (1. - 0.5 * TMath::Sqrt(K3));
  //  double K1 = K2*TMath::Sqrt(K3)/(2*TMath::ATan(TMath::Sqrt(K3)));
  double sqK3 = TMath::Sqrt(K3);
  double ATsqK3 = 0.5 / TMath::ATan(sqK3);
  double Y1 = lambda - w / h / 2;
  double Y2 = Y1 + w / h;
  double X1 = K2 * Y1;
  double X2 = K2 * Y2;
  double Z1 = sqK3 * TMath::TanH(X1);
  double Z2 = sqK3 * TMath::TanH(X2);
  double val = ATsqK3 * (TMath::ATan(Z2) - TMath::ATan(Z1));
  return val;
}


void  StTpcRSMaker::Print(Option_t* /* option */) const
{
  PrPP(Print, max_sectors);
  PrPP(Print, St_tpcPadConfigC::instance()->numberOfRows(1));
  PrPP(Print, St_tpcPadConfigC::instance()->numberOfRows(20));
  PrPP(Print, St_tpcPadConfigC::instance()->numberOfInnerRows(20));
  PrPP(Print, max_pads);
  PrPP(Print, St_TpcResponseSimulatorC::instance()->W());// = 26.2);//*eV
  PrPP(Print, St_TpcResponseSimulatorC::instance()->Cluster());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->longitudinalDiffusion());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->transverseDiffusion());
  //  PrPP(Print, Gain);
  PrPP(Print, max_timebins);
  PrPP(Print, numberOfInnerSectorAnodeWires);
  PrPP(Print, firstInnerSectorAnodeWire);
  PrPP(Print, lastInnerSectorAnodeWire);
  PrPP(Print, numberOfOuterSectorAnodeWires);
  PrPP(Print, firstOuterSectorAnodeWire);
  PrPP(Print, lastOuterSectorAnodeWire);
  PrPP(Print, anodeWirePitch);
  PrPP(Print, St_TpcResponseSimulatorC::instance()->OmegaTau()); // tan of Lorentz angle
  PrPP(Print, St_TpcResponseSimulatorC::instance()->NoElPerAdcI());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->NoElPerAdcO());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->NoElPerAdcX());
  PrPP(Print, anodeWireRadius);
  PrPP(Print, St_TpcResponseSimulatorC::instance()->AveragePedestal());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->AveragePedestalRMS());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->AveragePedestalRMSX());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->FanoFactor());

  for (int sector = 1; sector <= 24; sector++) {
    PrPP(Print, innerSectorAnodeVoltage[sector-1]);
    PrPP(Print, outerSectorAnodeVoltage[sector-1]);
  }

  PrPP(Print, St_TpcResponseSimulatorC::instance()->K3IP());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->K3IR());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->K3OP());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->K3OR());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->SigmaJitterTI());
  PrPP(Print, St_TpcResponseSimulatorC::instance()->SigmaJitterTO());
}


void StTpcRSMaker::DigitizeSector(int sector, tpcrs::DigiData& digi_data)
{
  static int iBreak = 0;
  static int AdcCut = 500;

  SignalSum_t* SignalSum = GetSignalSum(sector);
  double ped    = 0;
  int adc = 0;
  int index = 0;
  double gain = 1;
  int row, pad, bin;

  for (row = 1;  row <= St_tpcPadConfigC::instance()->numberOfRows(sector); row++) {
    int nPadsPerRow = St_tpcPadConfigC::instance()->padsPerRow(sector, row);
    double pedRMS = St_TpcResponseSimulatorC::instance()->AveragePedestalRMS();

    if (St_tpcAltroParamsC::instance()->N(sector - 1) > 0) {
      if (! (St_tpcPadConfigC::instance()->iTPC(sector) && St_tpcPadConfigC::instance()->IsRowInner(sector, row))) {
        pedRMS = St_TpcResponseSimulatorC::instance()->AveragePedestalRMSX();
      }
    }

#ifdef __DEBUG__
    float AdcSumBeforeAltro = 0, AdcSumAfterAltro = 0;
#endif

    for (pad = 1; pad <= nPadsPerRow; pad++) {
      gain = St_tpcPadGainT0BC::instance()->Gain(sector, row, pad);

      if (gain <= 0.0) continue;

      ped    = St_TpcResponseSimulatorC::instance()->AveragePedestal();
      static  short ADCs[max_timebins];
      static  short IDTs[max_timebins];
      memset(ADCs, 0, sizeof(ADCs));
      memset(IDTs, 0, sizeof(IDTs));
      int NoTB = 0;
      index = max_timebins * ((row - 1) * max_pads + pad - 1);

      for (bin = 0; bin < max_timebins; bin++, index++) {
        //	int index= max_timebins*((row-1)*max_pads+pad-1)+bin;
        // Digits : gain + ped
        //  GG TF1F *ff = new TF1F("ff","TMath::Sqrt(4.76658e+01*TMath::Exp(-2.87987e-01*(x-1.46222e+01)))",21,56)
        double pRMS = pedRMS;

        if (pRMS > 0) {
          adc = (int) (SignalSum[index].Sum / gain + gRandom->Gaus(ped, pRMS));
          adc = adc - (int) ped;
        }
        else
          adc = (int) (SignalSum[index].Sum / gain);

        if (adc > 1023) adc = 1023;

        if (adc < 1) continue;

        SignalSum[index].Adc = adc;
        NoTB++;
        ADCs[bin] = adc;
        IDTs[bin] = SignalSum[index].TrackId;
#ifdef __DEBUG__
        if (adc > 3 * pedRMS)	AdcSumBeforeAltro += adc;

        if (Debug() > 11 && SignalSum[index].Sum > 0) {
          LOG_INFO << "digi R/P/T/I = " << row << " /\t" << pad << " /\t" << bin << " /\t" << index
                   << "\tSum/Adc/TrackId = " << SignalSum[index].Sum << " /\t"
                   << SignalSum[index].Adc << " /\t" << SignalSum[index].TrackId << '\n';
        }
#endif
      }

      if (! NoTB) continue;

      if (St_tpcAltroParamsC::instance()->N(sector - 1) >= 0 && ! mAltro) {
        mAltro = new Altro(max_timebins, ADCs);

        if (St_tpcAltroParamsC::instance()->N(sector - 1) > 0) { // Tonko 06/25/08
          //      ConfigAltro(ONBaselineCorrection1, ONTailcancellation, ONBaselineCorrection2, ONClipping, ONZerosuppression)
          mAltro->ConfigAltro(                    0,                  1,                     0,          1,                 1);
          //       ConfigBaselineCorrection_1(int mode, int ValuePeDestal, int *PedestalMem, int polarity)
          //altro->ConfigBaselineCorrection_1(4, 0, PedestalMem, 0);  // Tonko 06/25/08
          mAltro->ConfigTailCancellationFilter(St_tpcAltroParamsC::instance()->K1(),
                                               St_tpcAltroParamsC::instance()->K2(),
                                               St_tpcAltroParamsC::instance()->K3(), // K1-3
                                               St_tpcAltroParamsC::instance()->L1(),
                                               St_tpcAltroParamsC::instance()->L2(),
                                               St_tpcAltroParamsC::instance()->L3());// L1-3
        }
        else {
          mAltro->ConfigAltro(0, 0, 0, 1, 1);
        }

        mAltro->ConfigZerosuppression(St_tpcAltroParamsC::instance()->Threshold(),
                                      St_tpcAltroParamsC::instance()->MinSamplesaboveThreshold(),
                                      0, 0);
        mAltro->PrintParameters();
      }

      if (mAltro) {
        mAltro->RunEmulation();
        NoTB = 0;
        int ADCsum = 0;

        for (int i = 0; i < max_timebins; i++) {
          if (ADCs[i] && ! mAltro->ADCkeep[i]) {ADCs[i] = 0;}

          if (ADCs[i]) {
            NoTB++;
            ADCsum += ADCs[i];
#ifdef __DEBUG__
            if (ADCs[i] > 3 * pedRMS) AdcSumAfterAltro += ADCs[i];

            if (Debug() > 12) {
              LOG_INFO << "Altro R/P/T/I = " << row << " /\t" << pad << " /\t" << i
                       << "\tAdc/TrackId = " << ADCs[i] << " /\t" << IDTs[i] << '\n';
            }
#endif
          }
          else {IDTs[i] = 0;}
        }
#ifdef __DEBUG__
        if (ADCsum > AdcCut) {
          iBreak++;
        }
#endif
      }
      else {
        if (St_tpcAltroParamsC::instance()->N(sector - 1) < 0) NoTB = AsicThresholds(ADCs);
      }

      if (NoTB > 0) {
        digi_data.Add(sector, row, pad, ADCs, IDTs, max_timebins);
      }
    } // pads

#ifdef __DEBUG__
    if (Debug() > 10) {
      LOG_INFO << "row = " << row << "\tAdcSumBeforeAltro = " << AdcSumBeforeAltro << "\tAdcSumAfterAltro = " << AdcSumAfterAltro << '\n';
    }
#endif
  } // row
}


int StTpcRSMaker::AsicThresholds(short* ADCs)
{
  int t1 = 0;
  int nSeqLo = 0;
  int nSeqHi = 0;
  int noTbleft = 0;

  for (unsigned int tb = 0; tb < max_timebins; tb++) {
    if (ADCs[tb] <= St_asic_thresholdsC::instance()->thresh_lo()) {
      if (! t1) ADCs[tb] = 0;
      else {
        if (nSeqLo <= St_asic_thresholdsC::instance()->n_seq_lo() ||
            nSeqHi <= St_asic_thresholdsC::instance()->n_seq_hi())
          for (unsigned int t = t1; t <= tb; t++) ADCs[t] = 0;
        else noTbleft += nSeqLo;
      }

      t1 = nSeqLo = nSeqHi = 0;
    }

    nSeqLo++;

    if (! t1) t1 = tb;

    if (ADCs[tb] > St_asic_thresholdsC::instance()->thresh_hi()) {nSeqHi++;}
  }

  return noTbleft;
}


double StTpcRSMaker::InducedCharge(double s, double h, double ra, double Va, double &t0)
{
  // Calculate variation of induced charge due to different arrived angles
  // alpha = -26 and -70 degrees
  LOG_INFO << "wire spacing = " << s << " cm"
           << "\tcathode anode gap = " << h << " cm"
           << "\tanode wire radius = " << ra << " cm"
           << "\tpotential on anode wire = " << Va << " V\n";
  const double B  = 30e-3; // 1/V
  const double E0 = 20e3; // V/cm
  const double mu = 2.26; // cm**2/V/sec CH4+ mobility
  // const double mu = 1.87; // cm**2/V/sec Ar+ mobility
  double alpha[2] = {-26., -70.};
  double pi = TMath::Pi();
  // E.Mathieson (3.2b), V.Chernyatin said that it should be used this (Weber ) approximation 07/09/08
  double rc = s / (2 * pi) * TMath::Exp(pi * h / s); LOG_INFO << "rc(Cylinder approx) = " << rc << " cm\n";
  //  double rc = 4*h/pi; LOG_INFO << "rc = " << rc << " cm\n";   // E.Mathieson (4.3), no valid for our case
  double C  = 1. / (2 * TMath::Log(rc / ra)); LOG_INFO << "C = " << C << '\n';
  double E  = 2 * pi * C * Va / s; LOG_INFO << "E = " << E << " V/cm\n";
  // Gain variation: M = M0*(1 - k*cos(2*alpha))
  double k = 2 * B / 3.*TMath::Power((pi / E0 / s), 2) * TMath::Power(C * Va, 3); LOG_INFO << "k = " << k << '\n';
  // Induced charge variation
  t0 = ra * ra / (4 * mu * C * Va);
  LOG_INFO << "t0 = " << 1e9 * t0 << " ns\n";                                   // E.Mathieson (2.10)
  double Tav = t0 * h / s / (2 * pi * C);  LOG_INFO << "Tav = " << 1e9 * Tav << " ns\n";
  //  double t = 5*55e-9;             LOG_INFO << "t = " << 1e9*t << " ns\n";
  double t = 180e-9;             LOG_INFO << "t = " << 1e9 * t << " ns\n";
  double rp = TMath::Sqrt(1. + t / t0); LOG_INFO << "r' = " << rp << '\n';
  // qc = rp*ra*sin(alpha)/(2*h) + C/2*log(1 + t/t0) = A*sin(alpha) + B
  double Aconstant = rp * ra / (2 * h);        LOG_INFO << "Aconstant = " << Aconstant << '\n';
  double Bconstant = C / 2 * TMath::Log(1 + t / t0); LOG_INFO << "Bconstant = " << Bconstant << '\n';
  double Gains[2];

  for (int i = 0; i < 2; i++) {
    Gains[i] = Aconstant * TMath::Sin(pi / 180 * alpha[i]) + Bconstant;
    LOG_INFO << "Gain = " << Gains[i] << " at alpha = " << alpha[i] << " degree\n";
  }

  double GainsAv = TMath::Sqrt(Gains[0] * Gains[1]);
  double r = 0;

  for (int i = 0; i < 2; i++) {
    r = TMath::Log(Gains[i] / GainsAv); LOG_INFO << "Relative gain " << r << " at alpha = " << alpha[i] << '\n';
  }

  return r;
}


int StTpcRSMaker::SearchT(const void* elem1, const void** elem2)
{
  g2t_tpc_hit_st* value1 = (g2t_tpc_hit_st*) elem1;
  g2t_tpc_hit_st* value2 = (g2t_tpc_hit_st*) *elem2;

  // sectors
  if ((value1->volume_id % 100000) / 100 != (value2->volume_id % 100000) / 100)
    return (value1->volume_id % 100000) / 100 - (value2->volume_id % 100000) / 100;

  // track id
  if (value1->track_p != value2->track_p)
    return value1->track_p - value2->track_p;

  // pad rows
  //  if (value1->volume_id%100 != value2->volume_id%100) return value1->volume_id%100 - value2->volume_id%100;
  // track length
  return (int) 100 * (value1->length - value2->length);
}


int StTpcRSMaker::CompareT(const void** elem1, const void** elem2)
{
  return SearchT(*elem1, elem2);
}


double StTpcRSMaker::fei(double t, double t0, double T)
{
  static const double xmaxt = 708.39641853226408;
  static const double xmaxD  = xmaxt - TMath::Log(xmaxt);
  double t01 = xmaxD, t11 = xmaxD;

  if (T > 0) {t11 = (t + t0) / T;}

  if (t11 > xmaxD) t11 = xmaxD;

  if (T > 0) {t01 = t0 / T;}

  if (t01 > xmaxD) t01  = xmaxD;

  return TMath::Exp(-t11) * (ROOT::Math::expint(t11) - ROOT::Math::expint(t01));
}


double StTpcRSMaker::shapeEI(double* x, double* par)  // does not work. It is needed to 1/s
{
  double t  = x[0];
  double value = 0;

  if (t <= 0) return value;

  double t0    = par[0];
  double T1 = par[1]; // tau_I
  double T2 = par[3]; // tau_C

  if (TMath::Abs((T1 - T2) / (T1 + T2)) < 1e-7) {
    return TMath::Max(0., (t + t0) / T1 * fei(t, t0, T1) + TMath::Exp(-t / T1) - 1);
  }

  if (T2 <= 0) return fei(t, t0, T1);

  if (T1 <= 0) return 0;

  return T1 / (T1 - T2) * (fei(t, t0, T1) - fei(t, t0, T2));
}


double StTpcRSMaker::shapeEI3(double* x, double* par)  // does not work. It is needed to 1/s
{
  double t  = x[0];
  double value = 0;

  if (t <= 0) return value;

  double t0    = par[0];
  double tau_F = par[1];
  double tau_P = par[2];
  double tau_I = par[3];
  double tau_C = par[5];
  double d =   1. / tau_P;
  double a[3] = {- 1. / tau_I, - 1. / tau_F, 0};
  double A[3] = {(a[0] + d) / (a[0] - a[1]), (a[1] + d) / (a[1] - a[0]), 0};
  int N = 2;

  if (tau_C > 0) {
    N = 3;
    a[2] = -1. / tau_C;
    A[0] = (a[0] + d) / a[0] / (a[0] - a[1]) / (a[0] - a[2]);
    A[1] = (a[1] + d) / a[1] / (a[1] - a[0]) / (a[1] - a[2]);
    A[2] = (a[2] + d) / a[2] / (a[2] - a[0]) / (a[2] - a[1]);
  }

  for (int i = 0; i < N; i++) {
    value += A[i] * TMath::Exp(a[i] * (t + t0)) * (ROOT::Math::expint(-a[i] * (t + t0)) - ROOT::Math::expint(-a[i] * t0));
  }

  return value;
}


double StTpcRSMaker::shapeEI_I(double* x, double* par)   //Integral of shape over time bin
{
  static double sqrt2 = TMath::Sqrt(2.);
  double TimeBinWidth = par[2];
  double norm = par[3];
  double t1 = TimeBinWidth * (x[0] - 0.5);
  double t2 = t1 + TimeBinWidth;
  int io = (int) par[4];
  assert(io >= 0 && io <= 1);
  return sqrt2 * fgTimeShape0[io]->Integral(t1, t2) / norm;
}


double StTpcRSMaker::shapeEI3_I(double* x, double* par)   //Integral of shape over time bin
{
  static double sqrt2 = TMath::Sqrt(2.);
  double TimeBinWidth = par[4];
  double norm = par[5];
  double t1 = TimeBinWidth * (x[0] - 0.5);
  double t2 = t1 + TimeBinWidth;
  int io = (int) par[6];
  assert(io >= 0 && io <= 1);
  return sqrt2 * fgTimeShape3[io]->Integral(t1, t2) / norm;
}


SignalSum_t*  StTpcRSMaker::GetSignalSum(int sector)
{
  if (! m_SignalSum)
    m_SignalSum = (SignalSum_t*) malloc(St_tpcPadConfigC::instance()->numberOfRows(sector) * max_pads * max_timebins * sizeof(SignalSum_t));

  return m_SignalSum;
}


SignalSum_t*  StTpcRSMaker::ResetSignalSum(int sector)
{
  GetSignalSum(sector);
  memset (m_SignalSum, 0, St_tpcPadConfigC::instance()->numberOfRows(sector)*max_pads * max_timebins * sizeof(SignalSum_t));
  return m_SignalSum;
}


double StTpcRSMaker::polya(double* x, double* par)
{
  return TMath::GammaDist(x[0], par[0], par[1], par[2]);
}


double StTpcRSMaker::Ec(double* x, double* p)
{
  if (x[0] < p[0] / 2 || x[0] > 3.064 * p[0]) return 0;

  if (x[0] < p[0]) return 1;

  return TMath::Power(p[0] / x[0], 4);
}


TF1* StTpcRSMaker::StTpcRSMaker::fEc(double w)
{
  TF1* f = new TF1("Ec", Ec, 0, 3.064 * w, 1);
  f->SetParameter(0, w);
  return f;
}


bool StTpcRSMaker::TrackSegment2Propagate(g2t_tpc_hit_st* tpc_hitC, const g2t_vertex_st* geantVertex, HitPoint_t &TrackSegmentHits)
{
  static int iBreak = 0;

  if (! tpc_hitC) {
    iBreak++;
  }

  if (tpc_hitC->de > 0) {
    mNSplittedHits = 0;
  }
  else if (! mNSplittedHits) {
    mNSplittedHits++;
  }

  int volId = tpc_hitC->volume_id % 10000;
  int sector = volId / 100;
  static StGlobalCoordinate coorG;    // ideal
  TrackSegmentHits.xyzG =
    StGlobalCoordinate(tpc_hitC->x[0], tpc_hitC->x[1], tpc_hitC->x[2]);  PrPP(Make, TrackSegmentHits.xyzG);
  coorG = TrackSegmentHits.xyzG;
  static StTpcLocalCoordinate  coorLT;  // before do distortions
  static StTpcLocalDirection  dirLT, BLT;
  // calculate row
  static StTpcLocalSectorCoordinate coorS;
  static StTpcCoordinateTransform transform(gStTpcDb);
  transform(coorG, coorS, sector, 0); PrPP(Make, coorS);
  int row = coorS.fromRow();
  transform(coorG, coorLT, sector, row); PrPP(Make, coorLT);
  int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;
  TrackSegmentHits.TrackId  = tpc_hitC->track_p;
  TrackSegmentHits.tpc_hitC = tpc_hitC;

  if (ClusterProfile) {
    checkList[io][0]->Fill(TrackSegmentHits.tpc_hitC->x[2], TMath::Abs(TrackSegmentHits.tpc_hitC->de));
    checkList[io][1]->Fill(TrackSegmentHits.tpc_hitC->x[2],           TrackSegmentHits.tpc_hitC->ds );
  }

  TrackSegmentHits.sMin = TrackSegmentHits.s - TrackSegmentHits.tpc_hitC->ds;
  TrackSegmentHits.sMax = TrackSegmentHits.s;

  if (TrackSegmentHits.sMin < msMin) msMin = TrackSegmentHits.sMin;

  if (TrackSegmentHits.sMax > msMax) msMax = TrackSegmentHits.sMax;

  // move up, calculate field at center of TPC
  static float BFieldG[3];
  StarMagField::Instance()->BField(tpc_hitC->x, BFieldG);
  // distortion and misalignment
  // replace pxy => direction and try linear extrapolation
  StThreeVector<double>       pxyzG(tpc_hitC->p[0], tpc_hitC->p[1], tpc_hitC->p[2]);
  StGlobalDirection    dirG(pxyzG.unit());                                    PrPP(Make, dirG);
  StGlobalDirection    BG(BFieldG[0], BFieldG[1], BFieldG[2]);                PrPP(Make, BG);
  transform( dirG,  dirLT, sector, row);                                      PrPP(Make, dirLT);
  transform(   BG,    BLT, sector, row);                                      PrPP(Make, BLT);

  // Distortions
  if (TESTBIT(m_Mode, kDistortion) && StMagUtilities::Instance()) {
    float pos[3] = {(float ) coorLT.position().x(), (float ) coorLT.position().y(), (float ) coorLT.position().z()};
    float posMoved[3];
    StMagUtilities::Instance()->DoDistortion(pos, posMoved, sector); // input pos[], returns posMoved[]
    StThreeVector<double> position(posMoved[0], posMoved[1], posMoved[2]);
    coorLT.setPosition(position);        // after do distortions
    transform(coorLT, TrackSegmentHits.xyzG);                PrPP(Make, coorLT);
  }

  // end of distortion
  transform(coorLT, TrackSegmentHits.coorLS); PrPP(Make, TrackSegmentHits.coorLS);
  transform( dirLT, TrackSegmentHits.dirLS);  PrPP(Make, TrackSegmentHits.dirLS);
  transform(   BLT, TrackSegmentHits.BLS);    PrPP(Make, TrackSegmentHits.BLS);

  double tof = geantVertex->ge_tof;
  //	if (! TESTBIT(m_Mode, kNoToflight))
  tof += tpc_hitC->tof;
  double driftLength = TrackSegmentHits.coorLS.position().z() + tof * gStTpcDb->DriftVelocity(sector); // ,row);

  if (driftLength > -1.0 && driftLength <= 0) {
    if ((row >  St_tpcPadConfigC::instance()->numberOfInnerRows(sector) && driftLength > -gStTpcDb->WirePlaneGeometry()->outerSectorAnodeWirePadPlaneSeparation()) ||
        (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector) && driftLength > -gStTpcDb->WirePlaneGeometry()->innerSectorAnodeWirePadPlaneSeparation()))
      driftLength = TMath::Abs(driftLength);
  }

  TrackSegmentHits.coorLS.position().setZ(driftLength); PrPP(Make, TrackSegmentHits.coorLS);
  // useT0, don't useTau
  transform(TrackSegmentHits.coorLS, TrackSegmentHits.Pad, false, false); // don't use T0, don't use Tau
  PrPP(Make, TrackSegmentHits.Pad);

  return true;
}


void StTpcRSMaker::GenerateSignal(HitPoint_t &TrackSegmentHits, int sector, int rowMin, int rowMax, double sigmaJitterT, double sigmaJitterX)
{
  static StTpcCoordinateTransform transform(gStTpcDb);
  SignalSum_t* SignalSum = GetSignalSum(sector);

  for (int row = rowMin; row <= rowMax; row++) {
    if (St_tpcPadConfigC::instance()->numberOfRows(sector) == 45) { // ! iTpx
      if ( ! StDetectorDbTpcRDOMasks::instance()->isRowOn(sector, row)) continue;

      if ( ! St_tpcAnodeHVavgC::instance()->livePadrow(sector, row))  continue;
    }

    int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;
    StTpcLocalSectorCoordinate xyzW(xOnWire, yOnWire, zOnWire, sector, row);
    static StTpcPadCoordinate Pad;
    transform(xyzW, Pad, false, false); // don't use T0, don't use Tau
    float bin = Pad.timeBucket();//L  - 1; // K
    int binT = tpcrs::irint(bin); //L bin;//K tpcrs::irint(bin);// J bin; // I tpcrs::irint(bin);

    if (binT < 0 || binT >= max_timebins) continue;

    double dT = bin -  binT + St_TpcResponseSimulatorC::instance()->T0offset();
    dT += (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ?
          St_TpcResponseSimulatorC::instance()->T0offsetI() :
          St_TpcResponseSimulatorC::instance()->T0offsetO();

    if (sigmaJitterT) dT += gRandom->Gaus(0, sigmaJitterT); // #1

    double dely      = {transform.yFromRow(sector, row) - yOnWire};
    double localYDirectionCoupling = mChargeFraction[io][sector - 1]->GetSaveL(&dely);

    if (ClusterProfile) {
      checkList[io][10]->Fill(TrackSegmentHits.xyzG.position().z(), localYDirectionCoupling);
    }

    if (localYDirectionCoupling < minSignal) continue;

    float padX = Pad.pad();
    int CentralPad = tpcrs::irint(padX);

    if (CentralPad < 1) continue;

    int PadsAtRow = St_tpcPadConfigC::instance()->numberOfPadsAtRow(sector, row);

    if (CentralPad > PadsAtRow) continue;

    int DeltaPad = tpcrs::irint(mPadResponseFunction[io][sector - 1]->GetXmax()) + 1;
    int padMin = TMath::Max(CentralPad - DeltaPad, 1);
    int padMax = std::min(CentralPad + DeltaPad, PadsAtRow);
    int Npads = std::min(padMax - padMin + 1, static_cast<int>(kPadMax));
    double xPadMin = padMin - padX;
    static double XDirectionCouplings[kPadMax];
    static double TimeCouplings[kTimeBacketMax];
    mPadResponseFunction[io][sector - 1]->GetSaveL(Npads, xPadMin, XDirectionCouplings);

    //	      double xPad = padMin - padX;
    for (int pad = padMin; pad <= padMax; pad++) {
      double gain = QAv * mGainLocal;
      double dt = dT;

      //		if (St_tpcPadConfigC::instance()->numberOfRows(sector) ==45 && ! TESTBIT(m_Mode, kGAINOAtALL)) {
      if (! TESTBIT(m_Mode, kGAINOAtALL)) {
        gain   *= St_tpcPadGainT0BC::instance()->Gain(sector, row, pad);

        if (gain <= 0.0) continue;

        dt -= St_tpcPadGainT0BC::instance()->T0(sector, row, pad);
      }

      if (ClusterProfile) {
        checkList[io][12]->Fill(TrackSegmentHits.xyzG.position().z(), gain);
        hist[4][1]->Fill(sector, row, gain);
      }

      //		double localXDirectionCoupling = localXDirectionCouplings[pad-padMin];
      double localXDirectionCoupling = gain * XDirectionCouplings[pad - padMin];

      if (localXDirectionCoupling < minSignal) continue;

      if (ClusterProfile) {
        checkList[io][13]->Fill(TrackSegmentHits.xyzG.position().z(), localXDirectionCoupling);
      }

      double XYcoupling = localYDirectionCoupling * localXDirectionCoupling;

      if (ClusterProfile) {
        checkList[io][14]->Fill(TrackSegmentHits.xyzG.position().z(), XYcoupling);
      }

      if (XYcoupling < minSignal)  continue;

      int bin_low  = TMath::Max(0, binT + tpcrs::irint(dt + mShaperResponse->GetXmin() - 0.5));
      int bin_high = std::min(max_timebins - 1, binT + tpcrs::irint(dt + mShaperResponse->GetXmax() + 0.5));
      int index = max_timebins * ((row - 1) * max_pads + pad - 1) + bin_low;
      int Ntbks = std::min(bin_high - bin_low + 1, static_cast<int>(kTimeBacketMax));
      double tt = -dt + (bin_low - binT);
      mShaperResponse->GetSaveL(Ntbks, tt, TimeCouplings);

      for (int itbin = bin_low; itbin <= bin_high; itbin++, index++) {
        double signal = XYcoupling * TimeCouplings[itbin - bin_low];

        if (signal < minSignal)  continue;

#ifdef __DEBUG__
        static int iBreak = 0;

        if (TMath::IsNaN(signal) || TMath::IsNaN(SignalSum[index].Sum)) {
          iBreak++;
        }
#endif
        TotalSignalInCluster += signal;
        SignalSum[index].Sum += signal;

        if (ClusterProfile) {
          if (pad >= pad0 && pad < pad0 + kPadMax &&
              itbin >= tbk0 &&  itbin < tbk0 + kTimeBacketMax) {
            padsdE[pad - pad0]   += signal;
            tbksdE[itbin - tbk0] += signal;
          }
        }

        rowsdE[row - 1]     += signal;
        rowsdEH[row - 1]    += signal;

        if ( TrackSegmentHits.TrackId ) {
          if (! SignalSum[index].TrackId ) SignalSum[index].TrackId = TrackSegmentHits.TrackId;
          else  // switch TrackId, works only for 2 tracks, more tracks ?
            if ( SignalSum[index].TrackId != TrackSegmentHits.TrackId && SignalSum[index].Sum < 2 * signal)
              SignalSum[index].TrackId = TrackSegmentHits.TrackId;
        }

#ifdef __DEBUG__
        if (Debug() > 13 && (SignalSum[index].Sum > 0 || ! TMath::Finite(SignalSum[index].Sum)) ) {
          LOG_INFO << "simu row = " << TrackSegmentHits.tpc_hitC->volume_id % 100 << "\tR/P/T/I = " << row << " /\t" << pad << " /\t" << itbin << " /\t" << index
                   << "\tSum/Adc/TrackId = " << SignalSum[index].Sum << " /\t"
                   << SignalSum[index].Adc << " /\t" << SignalSum[index].TrackId
                   << "\tsignal = " << signal
                   << "\trow Min/Max = " << rowMin << "/" << rowMax
                   << '\n';

          if (! TMath::Finite(SignalSum[index].Sum)) {
            LOG_INFO << "Not Finite\n";
          }
        }
#endif
      } // time
    } // pad limits
  } // row limits
}


double StTpcRSMaker::dEdxCorrection(HitPoint_t &TrackSegmentHits)
{
  double dEdxCor = 1;

  if (m_TpcdEdxCorrection) {
    //    dEdxCor = -1;
    double dStep =  TMath::Abs(TrackSegmentHits.tpc_hitC->ds);
    dEdxY2_t CdEdx;
    memset (&CdEdx, 0, sizeof(dEdxY2_t));
    CdEdx.DeltaZ = 5.2;
    CdEdx.QRatio = -2;
    CdEdx.QRatioA = -2.;
    CdEdx.QSumA = 0;
    CdEdx.sector = TrackSegmentHits.Pad.sector();
    CdEdx.row    = TrackSegmentHits.Pad.row();
    CdEdx.pad    = tpcrs::irint(TrackSegmentHits.Pad.pad());
    CdEdx.edge   = CdEdx.pad;

    if (CdEdx.edge > 0.5 * St_tpcPadConfigC::instance()->numberOfPadsAtRow(CdEdx.sector, CdEdx.row))
      CdEdx.edge += 1 - St_tpcPadConfigC::instance()->numberOfPadsAtRow(CdEdx.sector, CdEdx.row);

    CdEdx.F.dE     = 1;
    CdEdx.F.dx     = dStep;
    CdEdx.xyz[0] = TrackSegmentHits.coorLS.position().x();
    CdEdx.xyz[1] = TrackSegmentHits.coorLS.position().y();
    CdEdx.xyz[2] = TrackSegmentHits.coorLS.position().z();
    double probablePad = St_tpcPadConfigC::instance()->numberOfPadsAtRow(CdEdx.sector, CdEdx.row) / 2;
    double pitch = (CdEdx.row <= St_tpcPadConfigC::instance()->numberOfInnerRows(CdEdx.sector)) ?
                     St_tpcPadConfigC::instance()->innerSectorPadPitch(CdEdx.sector) :
                     St_tpcPadConfigC::instance()->outerSectorPadPitch(CdEdx.sector);
    double PhiMax = TMath::ATan2(probablePad * pitch, St_tpcPadConfigC::instance()->radialDistanceAtRow(CdEdx.sector, CdEdx.row));
    CdEdx.PhiR   = TMath::ATan2(CdEdx.xyz[0], CdEdx.xyz[1]) / PhiMax;
    CdEdx.xyzD[0] = TrackSegmentHits.dirLS.position().x();
    CdEdx.xyzD[1] = TrackSegmentHits.dirLS.position().y();
    CdEdx.xyzD[2] = TrackSegmentHits.dirLS.position().z();
    CdEdx.ZdriftDistance = CdEdx.xyzD[2];
    CdEdx.zG      = CdEdx.xyz[2];

    if (St_trigDetSumsC::instance())	CdEdx.Zdc     = St_trigDetSumsC::instance()->zdcX();

    CdEdx.ZdriftDistance = TrackSegmentHits.coorLS.position().z(); // drift length
    St_tpcGasC* tpcGas = m_TpcdEdxCorrection->tpcGas();

    if (tpcGas)
      CdEdx.ZdriftDistanceO2 = CdEdx.ZdriftDistance * tpcGas->Struct()->ppmOxygenIn;

    if (! m_TpcdEdxCorrection->dEdxCorrection(CdEdx)) {
      dEdxCor = CdEdx.F.dE;
    }
  }

  return dEdxCor;
}
#undef PrPP
