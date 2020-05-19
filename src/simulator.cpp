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

#include "simulator.h"

#include "Math/SpecFuncMathMore.h"
#include "TRandom.h"
#include "TFile.h"
#include "tcl.h"

#include "particles/StParticleTable.hh"
#include "particles/StParticleDefinition.hh"
#include "tpcrs/configurator.h"
#include "altro.h"
#include "bichsel.h"
#include "dedx_correction.h"
#include "logger.h"
#include "mag_utilities.h"
#include "math_funcs.h"
#include "struct_containers.h"
#include "track_helix.h"


#define __STOPPED_ELECTRONS__
#define __DEBUG__
#if defined(__DEBUG__)
#define PrPP(A,B) if (Debug()%10 > 2) {LOG_INFO << "StTpcRSMaker::" << (#A) << "\t" << (#B) << " = \t" << (B) << '\n';}
#else
#define PrPP(A,B)
#endif
#define LASERINO 170
#define CHASRINO 171

//                                    Inner        Outer
static       double t0IO[2]   = {1.20868e-9, 1.43615e-9}; // recalculated in InducedCharge

TF1F StTpcRSMaker::fgTimeShape3[2] = {
  TF1F("TimeShape3Inner;Time [s];Signal", StTpcRSMaker::shapeEI3, 0, 1, 7),
  TF1F("TimeShape3Outer;Time [s];Signal", StTpcRSMaker::shapeEI3, 0, 1, 7)
};

TF1F StTpcRSMaker::fgTimeShape0[2] = {
  TF1F("TimeShape0Inner;Time [s];Signal", StTpcRSMaker::shapeEI, 0, 1, 7),
  TF1F("TimeShape0Outer;Time [s];Signal", StTpcRSMaker::shapeEI, 0, 1, 7)
};

using dEdxCorr = StTpcdEdxCorrection::Corrections;
using tpcrs::Cfg;

StTpcRSMaker::StTpcRSMaker(double e_cutoff, const char* name):
  min_signal_(1e-4),
  electron_range_(0.0055), // Electron Range(.055mm)
  electron_range_energy_(3000), // eV
  electron_range_power_(1.78), // sigma =  electron_range_*(eEnery/electron_range_energy_)**electron_range_power_
  max_electron_energy_(e_cutoff),
  max_sectors_(Cfg<tpcDimensions>().numberOfSectors),
  max_pads_(182),
  max_timebins_(512),
  options_(0),
  mdNdx(nullptr),
  mdNdxL10(nullptr),
  mdNdEL10(nullptr),
  mShaperResponses{
    std::vector<TF1F>(max_sectors_, TF1F("ShaperFuncInner;Time [bin];Signal", StTpcRSMaker::shapeEI_I, 0, 1, 7)),
    std::vector<TF1F>(max_sectors_, TF1F("ShaperFuncOuter;Time [bin];Signal", StTpcRSMaker::shapeEI_I, 0, 1, 7))
  },
  mChargeFraction{
    std::vector<TF1F>(max_sectors_, TF1F("ChargeFractionInner;Distance [cm];Signal", StTpcRSMaker::PadResponseFunc, 0, 1, 6)),
    std::vector<TF1F>(max_sectors_, TF1F("ChargeFractionOuter;Distance [cm];Signal", StTpcRSMaker::PadResponseFunc, 0, 1, 6))
  },
  mPadResponseFunction{
    std::vector<TF1F>(max_sectors_, TF1F("PadResponseFunctionInner;Distance [pads];Signal", StTpcRSMaker::PadResponseFunc, 0, 1, 6)),
    std::vector<TF1F>(max_sectors_, TF1F("PadResponseFunctionOuter;Distance [pads];Signal", StTpcRSMaker::PadResponseFunc, 0, 1, 6))
  },
  mPolya{
    TF1F("PolyaInner;x = G/G_0;signal", polya, 0, 10, 3),
    TF1F("PolyaOuter;x = G/G_0;signal", polya, 0, 10, 3)
  },
  mHeed("Ec", StTpcRSMaker::Ec, 0, 3.064 * Cfg<TpcResponseSimulator>().W, 1),
  m_TpcdEdxCorrection(dEdxCorr::kAll & ~dEdxCorr::kAdcCorrection & ~dEdxCorr::kAdcCorrectionMDF & ~dEdxCorr::kdXCorrection, Debug())
{
  //  SETBIT(options_,kHEED);
  SETBIT(options_, kBICHSEL); // Default is Bichsel
  SETBIT(options_, kdEdxCorr);
  SETBIT(options_, kDistortion);

  if (TESTBIT(options_, kBICHSEL)) {
    LOG_INFO << "StTpcRSMaker:: use H.Bichsel model for dE/dx simulation\n";

    TFile inner(tpcrs::Configurator::Locate("dNdE_Bichsel.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(tpcrs::Configurator::Locate("dNdx_Bichsel.root").c_str());
    mdNdx = (TH1D*) outer.Get("dNdx"); assert(mdNdx);
    mdNdx->SetDirectory(0);
  }
  else if (TESTBIT(options_, kHEED)) {
    LOG_INFO << "StTpcRSMaker:: use Heed model for dE/dx simulation\n";

    TFile inner(tpcrs::Configurator::Locate("dNdx_Heed.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(tpcrs::Configurator::Locate("dNdx_Heed.root").c_str());
    mdNdxL10 = (TH1D*) outer.Get("dNdxL10"); assert(mdNdxL10);
    mdNdxL10->SetDirectory(0);
  }
  else {LOG_INFO << "StTpcRSMaker:: use GEANT321 model for dE/dx simulation\n";}

  if (TESTBIT(options_, kDistortion)) {
    LOG_INFO << "StTpcRSMaker:: use Tpc distortion correction\n";
  }

  double TimeBinWidth = 1. / Cfg<starClockOnl>().frequency;
  /*
  select firstInnerSectorAnodeWire,lastInnerSectorAnodeWire,numInnerSectorAnodeWires,firstOuterSectorAnodeWire,lastOuterSectorAnodeWire,numOuterSectorAnodeWires from  Geometry_tpc.tpcWirePlanes;
  +---------------------------+--------------------------+--------------------------+---------------------------+--------------------------+--------------------------+
  | firstInnerSectorAnodeWire | lastInnerSectorAnodeWire | numInnerSectorAnodeWires | firstOuterSectorAnodeWire | lastOuterSectorAnodeWire | numOuterSectorAnodeWires |
  +---------------------------+--------------------------+--------------------------+---------------------------+--------------------------+--------------------------+
  |             53.2000000000 |           120.8000000000 |                      170 |            122.7950000000 |           191.1950000000 |                      172 |
  +---------------------------+--------------------------+--------------------------+---------------------------+--------------------------+--------------------------+
   */
  numberOfInnerSectorAnodeWires  = Cfg<tpcWirePlanes>().numInnerSectorAnodeWires;
  firstInnerSectorAnodeWire      = Cfg<tpcWirePlanes>().firstInnerSectorAnodeWire;
  lastInnerSectorAnodeWire       = Cfg<tpcWirePlanes>().lastInnerSectorAnodeWire;
  numberOfOuterSectorAnodeWires  = Cfg<tpcWirePlanes>().numOuterSectorAnodeWires;
  firstOuterSectorAnodeWire      = Cfg<tpcWirePlanes>().firstOuterSectorAnodeWire;
  lastOuterSectorAnodeWire       = Cfg<tpcWirePlanes>().lastOuterSectorAnodeWire;
  anodeWirePitch                 = Cfg<tpcWirePlanes>().anodeWirePitch;
  anodeWireRadius                = Cfg<tpcWirePlanes>().anodeWireRadius;
  float BFieldG[3];
  float xyz[3] = {0, 0, 0};
  StarMagField::Instance().BField(xyz, BFieldG);
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
        if (sector > 1 && std::abs(innerSectorAnodeVoltage[sector - 1] - innerSectorAnodeVoltage[sector - 2]) < 1) {
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
        if (sector > 1 && std::abs(outerSectorAnodeVoltage[sector - 1] - outerSectorAnodeVoltage[sector - 2]) < 1) {
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
    FuncParams_t params3{
      {"t0",    t0IO[io]},
      {"tauF",  Cfg<TpcResponseSimulator>().tauF},
      {"tauP",  Cfg<TpcResponseSimulator>().tauP},
      {"tauI",  Cfg<TpcResponseSimulator>().tauIntegration},
      {"width", TimeBinWidth},
      {"tauC",  0},
      {"io",    io}
    };

    FuncParams_t params0{
      {"t0",    t0IO[io]},
      {"tauF",  0},
      {"tauP",  0},
      {"tauI",  io ? Cfg<TpcResponseSimulator>().tauXO : Cfg<TpcResponseSimulator>().tauXI},
      {"width", TimeBinWidth},
      {"tauC",  0},
      {"io",    io}
    };

    // old electronics, intergation + shaper alltogether
    for (int i = 0; i != params3.size(); ++i) {
      fgTimeShape3[io].SetParName(i, params3[i].first.c_str());
      fgTimeShape3[io].SetParameter(i, params3[i].second);
    }
    fgTimeShape3[io].SetRange(timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth);
    params3[5].second = fgTimeShape3[io].Integral(timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth);

    // new electronics only integration
    params0[5].second = io ? Cfg<TpcResponseSimulator>().tauCO : Cfg<TpcResponseSimulator>().tauCI;
    for (int i = 0; i != params0.size(); ++i) {
      fgTimeShape0[io].SetParName(i, params0[i].first.c_str());
      fgTimeShape0[io].SetParameter(i, params0[i].second);
    }
    fgTimeShape0[io].SetRange(timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth);
    params0[5].second = fgTimeShape0[io].Integral(0, timeBinMax * TimeBinWidth);

    for (int sector = 1; sector <= max_sectors_; sector++) {
      //                            w       h        s       a       l  i
      //  double paramsI[6] = {0.2850, 0.2000,  0.4000, 0.0010, 1.1500, 0};
      //  double paramsO[6] = {0.6200, 0.4000,  0.4000, 0.0010, 1.1500, 0};
      double params[6]{
        io == kInner ? St_tpcPadConfigC::instance()->innerSectorPadWidth(sector) :  // w = width of pad
                       St_tpcPadConfigC::instance()->outerSectorPadWidth(sector),
        io == kInner ? Cfg<tpcWirePlanes>().innerSectorAnodeWirePadSep :            // h = Anode-Cathode gap
                       Cfg<tpcWirePlanes>().outerSectorAnodeWirePadSep,
        Cfg<tpcWirePlanes>().anodeWirePitch,                                        // s = wire spacing
        io == kInner ? Cfg<TpcResponseSimulator>().K3IP :
                       Cfg<TpcResponseSimulator>().K3OP,
        0,
        io == kInner ? St_tpcPadConfigC::instance()->innerSectorPadPitch(sector) :
                       St_tpcPadConfigC::instance()->outerSectorPadPitch(sector)
      };

      mPadResponseFunction[io][sector - 1].SetParameters(params);
      mPadResponseFunction[io][sector - 1].SetParNames("PadWidth", "Anode-Cathode gap", "wire spacing", "K3OP", "CrossTalk", "PadPitch");
      mPadResponseFunction[io][sector - 1].SetRange(-2.5, 2.5); // Cut tails
      mPadResponseFunction[io][sector - 1].Save(-4.5, 4.5, 0, 0, 0, 0);

      params[0] = io == kInner ? St_tpcPadConfigC::instance()->innerSectorPadLength(sector) :
                                 St_tpcPadConfigC::instance()->outerSectorPadLength(sector);
      params[3] = io == kInner ? Cfg<TpcResponseSimulator>().K3IR :
                                 Cfg<TpcResponseSimulator>().K3OR;
      params[5] = 1.;

      // Cut the tails
      double x_range = 2.5;
      for (; x_range > 1.5; x_range -= 0.05) {
        double r = mChargeFraction[io][sector - 1].Eval(x_range) / mChargeFraction[io][sector - 1].Eval(0);
        if (r > 1e-2) break;
      }

      mChargeFraction[io][sector - 1].SetParameters(params);
      mChargeFraction[io][sector - 1].SetParNames("PadLength", "Anode-Cathode gap", "wire spacing", "K3IR", "CrossTalk", "RowPitch");
      mChargeFraction[io][sector - 1].SetRange(-x_range, x_range);
      mChargeFraction[io][sector - 1].Save(-2.5, 2.5, 0, 0, 0, 0);

      //  TF1F *func = new TF1F("funcP","x*sqrt(x)/exp(2.5*x)",0,10);
      // see http://www4.rcf.bnl.gov/~lebedev/tec/polya.html
      // Gain fluctuation in proportional counters follows Polya distribution.
      // x = G/G_0
      // P(m) = m(m(x)**(m-1)*exp(-m*x)/Gamma(m);
      // original Polya  m = 1.5 (R.Bellazzini and M.A.Spezziga, INFN PI/AE-94/02).
      // Valeri Cherniatin (cherniat@bnlarm.bnl.gov) recomends m=1.38
      // Trs uses x**1.5/exp(x)
      // tss used x**0.5/exp(1.5*x)
      if (Cfg<tpcAltroParams>(sector - 1).N < 0) { // old TPC
        InitShaperFuncs(io, sector, mShaperResponses, StTpcRSMaker::shapeEI3_I, params3, timeBinMin, timeBinMax);
      } else {//Altro
        InitShaperFuncs(io, sector, mShaperResponses, StTpcRSMaker::shapeEI_I,  params0, timeBinMin, timeBinMax);
      }
    }
  }

  if (Debug()) Print();

  //  mPolya = new TF1F("Polya;x = G/G_0;signal","sqrt(x)/exp(1.5*x)",0,10); // original Polya
  //  mPolya = new TF1F("Polya;x = G/G_0;signal","pow(x,0.38)*exp(-1.38*x)",0,10); //  Valeri Cherniatin
  //   mPoly = new TH1D("Poly","polyaAvalanche",100,0,10);
  //if (gamma <= 0) gamma = 1.38;
  double gamma_inn = Cfg<TpcResponseSimulator>().PolyaInner;
  double gamma_out = Cfg<TpcResponseSimulator>().PolyaOuter;
  mPolya[kInner].SetParameters(gamma_inn, 0., 1. / gamma_inn);
  mPolya[kOuter].SetParameters(gamma_out, 0., 1. / gamma_out);

  // HEED function to generate Ec, default w = 26.2
  mHeed.SetParameter(0, Cfg<TpcResponseSimulator>().W);
}


StTpcRSMaker::~StTpcRSMaker()
{
  delete mdNdx;
  delete mdNdxL10;
  delete mdNdEL10;
}


void StTpcRSMaker::InitShaperFuncs(int io, int sector, std::array<std::vector<TF1F>, 2>& funcs,
  double (*shape)(double*, double*), FuncParams_t params, double timeBinMin, double timeBinMax)
{
  funcs[io][sector - 1].SetFunction(shape);
  funcs[io][sector - 1].SetRange(timeBinMin, timeBinMax);

  for (int i = 0; i != params.size(); ++i) {
    funcs[io][sector - 1].SetParName(i, params[i].first.c_str());
    funcs[io][sector - 1].SetParameter(i, params[i].second);
  }

  // Cut tails
  double t = timeBinMax;
  double ymax = funcs[io][sector - 1].Eval(0.5);

  for (; t > 5; t -= 1) {
    double r = funcs[io][sector - 1].Eval(t) / ymax;
    if (r > 1e-2) break;
  }

  funcs[io][sector - 1].SetRange(timeBinMin, t);
  funcs[io][sector - 1].Save(timeBinMin, t, 0, 0, 0, 0);
}


void StTpcRSMaker::Make(std::vector<tpcrs::GeantHit>& geant_hits, tpcrs::DigiData& digi_data)
{
  static int nCalls = 0;
  gRandom->SetSeed(2345 + nCalls++);

  double vminI = St_tpcGainCorrectionC::instance()->Struct(1)->min;
  double vminO = St_tpcGainCorrectionC::instance()->Struct(0)->min;

  // TODO: Confirm proper handling of empty input containers

  St_tpcGainCorrectionC::instance()->Struct(0)->min = -500;
  St_tpcGainCorrectionC::instance()->Struct(1)->min = -500;

  if (Debug()) {
    LOG_INFO << "Reset min for gain Correction to I/O\t"
             << St_tpcGainCorrectionC::instance()->Struct(1)->min
             << "\t"
             << St_tpcGainCorrectionC::instance()->Struct(0)->min
             << " (V)\n";
  }

  // sort
  // initialize original index locations
  int n_hits = geant_hits.size();
  std::vector<size_t> sorted_index(n_hits);
  std::iota(sorted_index.begin(), sorted_index.end(), 0);
  std::stable_sort(sorted_index.begin(), sorted_index.end(), [&geant_hits](size_t i1, size_t i2) {return geant_hits[i1] < geant_hits[i2];});

  int sortedIndex = 0;

  for (int sector = 1; sector <= max_sectors_; sector++) {
    int nHitsInTheSector = 0;
    std::vector<SignalSum_t> binned_charge(St_tpcPadConfigC::instance()->numberOfRows(sector) * max_pads_ * max_timebins_);

    // it is assumed that hit are ordered by sector, trackId, pad rows, and track length
    for (; sortedIndex < n_hits; sortedIndex++) {
      int indx = sorted_index[sortedIndex];

      tpcrs::GeantHit& geant_hit = geant_hits[indx];
      int volId = geant_hit.volume_id % 10000;
      int iSector = volId / 100;

      if (iSector != sector) {
        if (iSector < sector) {
          LOG_ERROR << "StTpcRSMaker::Make: geant_hits table has not been ordered by sector no. " << sector << '\n';
          assert( iSector > sector );
        }

        break;
      }

      if (geant_hit.volume_id <= 0 || geant_hit.volume_id > 1000000) continue;

      int parent_track_idx  = geant_hit.track_id;
      double mass = 0;

      int ipart      = geant_hit.particle_id;
      int charge     = 0;

      StParticleDefinition* particle = StParticleTable::instance()->findParticleByGeantId(ipart);

      if (particle) {
        mass = particle->mass();
        charge = particle->charge();
      }

      if (ipart == LASERINO || ipart == CHASRINO) {
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
      }
      // special treatment for electron/positron
      if (ipart == 2) charge =  101;
      if (ipart == 3) charge = -101;

      // Track segment to propagate
      double smin = 9999;
      double smax = -9999;
      int sIndex = sortedIndex;

      std::vector<HitPoint_t> TrackSegmentHits;
      BuildTrackSegments(sector, sorted_index, sortedIndex, geant_hits, TrackSegmentHits, smin, smax, sIndex);
      int nSegHits = TrackSegmentHits.size();

      if (!nSegHits) continue;

      if (Debug() >= 10) {
        PrPP(Make, nSegHits);

        for (int s = 0; s < nSegHits; s++) {
          LOG_INFO << "Seg[" << Form("%2i", s) << "]\tId " << TrackSegmentHits[s].TrackId << "\ts = " << TrackSegmentHits[s].s
               << "\tvolumeID :" <<  Form("%6i", TrackSegmentHits[s].tpc_hitC->volume_id) << "\t" << TrackSegmentHits[s].Pad
               << "\ts1/s2 = " << TrackSegmentHits[s].tpc_hitC->len - TrackSegmentHits[s].tpc_hitC->ds / 2
               << "\t" << TrackSegmentHits[s].tpc_hitC->len + TrackSegmentHits[s].tpc_hitC->ds / 2 << "\tds = " << TrackSegmentHits[s].tpc_hitC->ds
               << '\n';
        }
      }

      sortedIndex = sIndex - 1; // Irakli 05/06/19, reduce extra step in for loop

      double s = smin;
      memset (rowsdE, 0, sizeof(rowsdE));

      for (int iSegHits = 0; iSegHits < nSegHits && s < smax; iSegHits++) {
        tpcrs::GeantHit* tpc_hitC = TrackSegmentHits[iSegHits].tpc_hitC;
        tpc_hitC->digi.adc = 0;
        int row = TrackSegmentHits[iSegHits].coorLS.row;
        int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;
        // switch between Inner / Outer Sector paramters
        // Extra correction for simulation with respect to data
        int iowe = 0;
        if (sector > 12) iowe += 4;
        if (io)          iowe += 2;

        const float* AdditionalMcCorrection = Cfg<TpcResponseSimulator>().SecRowCorIW;
        const float* AddSigmaMcCorrection   = Cfg<TpcResponseSimulator>().SecRowSigIW;
        // Generate signal
        double sigmaJitterT = Cfg<TpcResponseSimulator>().SigmaJitterTI;
        double sigmaJitterX = Cfg<TpcResponseSimulator>().SigmaJitterXI;

        if (io) { // Outer
          sigmaJitterT = Cfg<TpcResponseSimulator>().SigmaJitterTO;
          sigmaJitterX = Cfg<TpcResponseSimulator>().SigmaJitterXO;
        }

        // Generate signal
        double Gain = GainCorrection(sector, row);

        double GainXCorrectionL = AdditionalMcCorrection[iowe] + row * AdditionalMcCorrection[iowe + 1];
        Gain *= std::exp(-GainXCorrectionL);
        double GainXSigma = AddSigmaMcCorrection[iowe] + row * AddSigmaMcCorrection[iowe + 1];

        if (GainXSigma > 0) Gain *= std::exp(gRandom->Gaus(0., GainXSigma));

        // dE/dx correction
        double dEdxCor = dEdxCorrection(TrackSegmentHits[iSegHits]);
        if (dEdxCor <= 0.) continue;

        dEdxCor *= GatingGridTransparency(TrackSegmentHits[iSegHits].Pad.timeBucket);

        if (dEdxCor < min_signal_) continue;

        // Initialize propagation
        // Magnetic field BField must be in kilogauss
        // kilogauss = 1e-1*tesla = 1e-1*(volt*second/meter2) = 1e-1*(1e-6*1e-3*1/1e4) = 1e-14
        TrackHelix track(TrackSegmentHits[iSegHits].dirLS.position,
                         TrackSegmentHits[iSegHits].coorLS.position,
                         TrackSegmentHits[iSegHits].BLS.position.z * 1e-14 * charge, 1);
#ifdef __DEBUG__
        if (Debug() > 11) PrPP(Make, track);
#endif
        double dStep =  std::abs(tpc_hitC->ds);
        double s_low   = -dStep / 2;
        double s_upper = s_low + dStep;
        double newPosition = s_low;
        // Propagate track to the pad row plane defined by the normal in this sector coordinate system
        static CoordTransform transform;
        Coords rowPlane{0, transform.yFromRow(TrackSegmentHits[iSegHits].Pad.sector, TrackSegmentHits[iSegHits].Pad.row), 0};
        double sR = track.pathLength(rowPlane, {0, 1, 0});

        if (sR < 1e10) {
          PrPP(Maker, sR);
          PrPP(Make, TrackSegmentHits[iSegHits].coorLS);
          TrackSegmentHits[iSegHits].coorLS.position = {track.at(sR).x, track.at(sR).y, track.at(sR).z};
          PrPP(Make, TrackSegmentHits[iSegHits].coorLS);
          PrPP(Make, TrackSegmentHits[iSegHits].Pad);
          transform.local_sector_to_hardware(TrackSegmentHits[iSegHits].coorLS, TrackSegmentHits[iSegHits].Pad, false, false); // don't use T0, don't use Tau
          PrPP(Make, TrackSegmentHits[iSegHits].Pad);
        }

        static const double m_e = .51099907e-3;
        static const double eV = 1e-9; // electronvolt in GeV
        double total_signal = 0;
        double gamma = std::pow(10., tpc_hitC->lgam) + 1;
        double betaGamma = std::sqrt(gamma * gamma - 1.);
        double bg = 0;
        double eKin = -1;
        Coords pxyzG{tpc_hitC->p[0], tpc_hitC->p[1], tpc_hitC->p[2]};

#ifdef __STOPPED_ELECTRONS__
        if (mass > 0) {
          bg = pxyzG.mag() / mass;

          // special case stopped electrons
          if (tpc_hitC->ds < 0.0050 && tpc_hitC->de < 0) {
            int Id    = tpc_hitC->track_id;
            int ipart = tpc_hitC->particle_id;

            if (ipart == 3) {
              eKin = -tpc_hitC->de;
              gamma = eKin / m_e + 1;
              bg = std::sqrt(gamma * gamma - 1.);
            }
          }
        }
#else
        if (mass > 0) bg = pxyzG.mag() / mass;
#endif
        if (bg > betaGamma) betaGamma = bg;

        gamma = std::sqrt(betaGamma * betaGamma + 1.);
        double Tmax;

        if (mass < 2 * m_e) {
          if (charge > 0) Tmax =       m_e * (gamma - 1);
          else            Tmax = 0.5 * m_e * (gamma - 1);
        }
        else {
          double r = m_e / mass;
          Tmax = 2 * m_e * betaGamma * betaGamma / (1 + 2 * gamma * r + r * r);
        }

        if (Tmax > max_electron_energy_) Tmax = max_electron_energy_;

        double OmegaTau = Cfg<TpcResponseSimulator>().OmegaTau *
                            TrackSegmentHits[iSegHits].BLS.position.z / 5.0; // from diffusion 586 um / 106 um at B = 0/ 5kG


        double driftLength = std::abs(TrackSegmentHits[iSegHits].coorLS.position.z);
        double D = 1. + OmegaTau * OmegaTau;
        double SigmaT = Cfg<TpcResponseSimulator>().transverseDiffusion * std::sqrt(driftLength / D);

        //	double SigmaL = Cfg<TpcResponseSimulator>().longitudinalDiffusion*std::sqrt(2*driftLength  );
        if (sigmaJitterX > 0) {SigmaT = std::sqrt(SigmaT * SigmaT + sigmaJitterX * sigmaJitterX);}

        double SigmaL     = Cfg<TpcResponseSimulator>().longitudinalDiffusion * std::sqrt(driftLength);
        double NoElPerAdc = Cfg<TpcResponseSimulator>().NoElPerAdc;

        if (NoElPerAdc <= 0) {
          if (St_tpcPadConfigC::instance()->iTPC(sector) && St_tpcPadConfigC::instance()->IsRowInner(sector, row)) {
            NoElPerAdc = Cfg<TpcResponseSimulator>().NoElPerAdcX; // iTPC
          }
          else if (St_tpcPadConfigC::instance()->IsRowInner(sector, row)) {
            NoElPerAdc = Cfg<TpcResponseSimulator>().NoElPerAdcI; // inner TPX
          }
          else {
            NoElPerAdc = Cfg<TpcResponseSimulator>().NoElPerAdcO; // outer TPX
          }
        }

#ifndef __NO_1STROWCORRECTION__
        if (row == 1) dEdxCor *= std::exp(Cfg<TpcResponseSimulator>().FirstRowC);
#endif
        double gain_local = Gain / dEdxCor / NoElPerAdc; // Account dE/dx calibration
        // end of dE/dx correction
        // generate electrons: No. of primary clusters per cm

        double NP;
        if (mdNdx || mdNdxL10) {
          NP = GetNoPrimaryClusters(betaGamma, charge); // per cm
#ifdef __DEBUG__
          if (NP <= 0.0) {
            continue;
          }
#endif
        }

        int nP = 0;
        double dESum = 0;
        double dSSum = 0;
        int   nTotal = 0;
        float dEr = 0;

        do {// Clusters
          float dS = 0;
          float dE = 0;
          static double cLog10 = std::log(10.);

          if (eKin >= 0.0) {
            if (eKin == 0.0) break;

            gamma = eKin / m_e + 1;
            bg = std::sqrt(gamma * gamma - 1.);
            Tmax = 0.5 * m_e * (gamma - 1);

            if (Tmax <= Cfg<TpcResponseSimulator>().W / 2 * eV) break;

            NP = GetNoPrimaryClusters(betaGamma, charge);
            dE = std::exp(cLog10 * mdNdEL10->GetRandom());
          }
          else {
            if (charge) {
              dS = - std::log(gRandom->Rndm()) / NP;

              if (mdNdEL10) dE = std::exp(cLog10 * mdNdEL10->GetRandom());
              else          dE = Cfg<TpcResponseSimulator>().W *
                                 gRandom->Poisson(Cfg<TpcResponseSimulator>().Cluster);
            }
            else { // charge == 0 geantino
              // for LASERINO assume dE/dx = 25 keV/cm;
              dE = 10; // eV
              dS = dE * eV / (std::abs(tpc_hitC->de / tpc_hitC->ds));
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

          if (dE < Cfg<TpcResponseSimulator>().W / 2 || E > Tmax) continue;

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
          if (dE > electron_range_energy_)
            xRange = electron_range_ * std::pow((dE + dEr) / electron_range_energy_, electron_range_power_);

          std::vector<float> rs = NumberOfElectronsInCluster(mHeed, dE, dEr);

          if (!rs.size()) continue;

          Coords xyzC = track.at(newPosition);

          double total_signal_in_cluster =
            LoopOverElectronsInCluster(rs, TrackSegmentHits[iSegHits], binned_charge, sector, row, xRange, xyzC, gain_local, SigmaT, SigmaL, OmegaTau);
          nTotal = rs.size();

          total_signal += total_signal_in_cluster;
        }
        while (true);   // Clusters

#ifdef __DEBUG__
        if (Debug() > 12) {
          LOG_INFO << "sIndex = " << sIndex << " volId = " << volId
                   << " dESum = " << dESum << " /\tdSSum " << dSSum << " /\t total_signal " << total_signal << '\n';
        }
#endif
        tpc_hitC->digi.de = dESum * eV;
        tpc_hitC->digi.ds = dSSum;
        tpc_hitC->digi.np = nP;
        tpc_hitC->digi.pad = TrackSegmentHits[iSegHits].Pad.pad;
        tpc_hitC->digi.timebin = TrackSegmentHits[iSegHits].Pad.timeBucket;

        nHitsInTheSector++;
      } // end do loop over segments for a given particle

      for (int iSegHits = 0; iSegHits < nSegHits; iSegHits++) {
        tpcrs::GeantHit* tpc_hitC = TrackSegmentHits[iSegHits].tpc_hitC;

        if (tpc_hitC->volume_id > 10000) continue;

        int row = tpc_hitC->volume_id % 100;
        tpc_hitC->digi.adc += rowsdE[row - 1];
      }
    }  // hits in the sector

    if (nHitsInTheSector) {
      DigitizeSector(sector, digi_data, binned_charge);

      if (Debug()) LOG_INFO << "StTpcRSMaker: Done with sector\t" << sector << " total no. of hit = " << nHitsInTheSector << '\n';
    }
  } // sector

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


void StTpcRSMaker::BuildTrackSegments(int sector, const std::vector<size_t>& sorted_index, int sortedIndex,
  std::vector<tpcrs::GeantHit>& geant_hits,
  std::vector<HitPoint_t>& segments, double& smin, double& smax, int& sIndex)
{
  int n_hits = sorted_index.size();

  if (Debug() > 13) LOG_INFO << "sortedIndex = " << sortedIndex << "\tn_hits = " << n_hits << '\n';

  segments.resize(100);
  HitPoint_t prev_segment;
  int parent_track_idx = 0;
  int TrackDirection = 0; // 0 - increase no of row, 1 - decrease no of. row.
  int num_segments = 0;

  for (sIndex = sortedIndex; sIndex < n_hits && num_segments < 100; sIndex++)
  {
    int indx = sorted_index[sIndex];
    tpcrs::GeantHit& geant_hit = geant_hits[indx];

    if ((geant_hit.volume_id % 10000) / 100 != sector) break;

    if (parent_track_idx > 0 && parent_track_idx != geant_hit.track_id) break;

    parent_track_idx = geant_hit.track_id;

    if (num_segments == 1) { // No Loopers !
      if (prev_segment.tpc_hitC->volume_id % 100 <= geant_hit.volume_id % 100) {
        TrackDirection = 0;
      }
      else {
        TrackDirection = 1;
      }
    }
    else if (num_segments > 1) {
      if ((! TrackDirection && prev_segment.tpc_hitC->volume_id % 100 > geant_hit.volume_id % 100) ||
          (  TrackDirection && prev_segment.tpc_hitC->volume_id % 100 < geant_hit.volume_id % 100))
        break;
    }

    if (Debug() > 13) LOG_INFO << "sIndex = " << sIndex << "\tindx = " << indx << "\ttpc_hitC = " << &geant_hit << '\n';

    HitPoint_t curr_segment;
    curr_segment.s = geant_hit.len;

    if (geant_hit.len == 0 && num_segments > 1) {
      curr_segment.s = prev_segment.s + curr_segment.tpc_hitC->ds;
    }

    TrackSegment2Propagate(geant_hit, curr_segment, smin, smax);

    if (curr_segment.Pad.timeBucket < 0 || curr_segment.Pad.timeBucket > max_timebins_) continue;

    segments[num_segments++] = curr_segment;
    prev_segment = curr_segment;
  }

  segments.resize(num_segments);
}


double StTpcRSMaker::GetNoPrimaryClusters(double betaGamma, int charge)
{
  double beta = betaGamma / std::sqrt(1.0 + betaGamma * betaGamma);
  double dNdx = 0;

  if      (mdNdx   ) dNdx = mdNdx->Interpolate(betaGamma);
  else if (mdNdxL10) dNdx = mdNdxL10->Interpolate(std::log10(betaGamma));

  double Q_eff = std::abs(charge % 100);

  if (Q_eff > 1)   {
    // Effective charge from GEANT ghion.F
    double w1 = 1.034 - 0.1777 * std::exp(-0.08114 * Q_eff);
    double w2 = beta * std::pow(Q_eff, -2. / 3.);
    double w3 = 121.4139 * w2 + 0.0378 * std::sin(190.7165 * w2);
    Q_eff    *= 1. - w1 * std::exp(-w3);
  }

  return Q_eff * Q_eff * dNdx;
}


double StTpcRSMaker::PadResponseFunc(double* x, double* par)
{
  double cross_talk = 0;
  double pad_response = 0;
  double X = par[5] * x[0];

  if (cross_talk > 0) {
    for (int i = -1; i <= 1; i++) {
      double xx = X + par[5] * i;

      if (i == 0) pad_response += (1. - 2.*cross_talk) * Gatti(&xx, par);
      else        pad_response +=          cross_talk  * Gatti(&xx, par);
    }
  }
  else pad_response = Gatti(&X, par);

  return pad_response;
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
  double K2 = M_PI_2 * (1. - 0.5 * std::sqrt(K3));
  //  double K1 = K2*std::sqrt(K3)/(2*std::atan(std::sqrt(K3)));
  double sqK3 = std::sqrt(K3);
  double ATsqK3 = 0.5 / std::atan(sqK3);
  double Y1 = lambda - w / h / 2;
  double Y2 = Y1 + w / h;
  double X1 = K2 * Y1;
  double X2 = K2 * Y2;
  double Z1 = sqK3 * std::tanh(X1);
  double Z2 = sqK3 * std::tanh(X2);
  double val = ATsqK3 * (std::atan(Z2) - std::atan(Z1));
  return val;
}


void  StTpcRSMaker::Print(Option_t* /* option */) const
{
  PrPP(Print, max_sectors_);
  PrPP(Print, St_tpcPadConfigC::instance()->numberOfRows(1));
  PrPP(Print, St_tpcPadConfigC::instance()->numberOfRows(20));
  PrPP(Print, St_tpcPadConfigC::instance()->numberOfInnerRows(20));
  PrPP(Print, max_pads_);
  PrPP(Print, Cfg<TpcResponseSimulator>().W);// = 26.2);//*eV
  PrPP(Print, Cfg<TpcResponseSimulator>().Cluster);
  PrPP(Print, Cfg<TpcResponseSimulator>().longitudinalDiffusion);
  PrPP(Print, Cfg<TpcResponseSimulator>().transverseDiffusion);
  //  PrPP(Print, Gain);
  PrPP(Print, max_timebins_);
  PrPP(Print, numberOfInnerSectorAnodeWires);
  PrPP(Print, firstInnerSectorAnodeWire);
  PrPP(Print, lastInnerSectorAnodeWire);
  PrPP(Print, numberOfOuterSectorAnodeWires);
  PrPP(Print, firstOuterSectorAnodeWire);
  PrPP(Print, lastOuterSectorAnodeWire);
  PrPP(Print, anodeWirePitch);
  PrPP(Print, Cfg<TpcResponseSimulator>().OmegaTau); // tan of Lorentz angle
  PrPP(Print, Cfg<TpcResponseSimulator>().NoElPerAdcI);
  PrPP(Print, Cfg<TpcResponseSimulator>().NoElPerAdcO);
  PrPP(Print, Cfg<TpcResponseSimulator>().NoElPerAdcX);
  PrPP(Print, anodeWireRadius);
  PrPP(Print, Cfg<TpcResponseSimulator>().AveragePedestal);
  PrPP(Print, Cfg<TpcResponseSimulator>().AveragePedestalRMS);
  PrPP(Print, Cfg<TpcResponseSimulator>().AveragePedestalRMSX);
  PrPP(Print, Cfg<TpcResponseSimulator>().FanoFactor);

  for (int sector = 1; sector <= 24; sector++) {
    PrPP(Print, innerSectorAnodeVoltage[sector-1]);
    PrPP(Print, outerSectorAnodeVoltage[sector-1]);
  }

  PrPP(Print, Cfg<TpcResponseSimulator>().K3IP);
  PrPP(Print, Cfg<TpcResponseSimulator>().K3IR);
  PrPP(Print, Cfg<TpcResponseSimulator>().K3OP);
  PrPP(Print, Cfg<TpcResponseSimulator>().K3OR);
  PrPP(Print, Cfg<TpcResponseSimulator>().SigmaJitterTI);
  PrPP(Print, Cfg<TpcResponseSimulator>().SigmaJitterTO);
}


void StTpcRSMaker::DigitizeSector(int sector, tpcrs::DigiData& digi_data, std::vector<SignalSum_t>& binned_charge)
{
  for (int row = 1;  row <= St_tpcPadConfigC::instance()->numberOfRows(sector); row++) {
    int nPadsPerRow = St_tpcPadConfigC::instance()->padsPerRow(sector, row);
    double pedRMS = Cfg<TpcResponseSimulator>().AveragePedestalRMS;

    if (Cfg<tpcAltroParams>(sector - 1).N > 0) {
      if (! (St_tpcPadConfigC::instance()->iTPC(sector) && St_tpcPadConfigC::instance()->IsRowInner(sector, row))) {
        pedRMS = Cfg<TpcResponseSimulator>().AveragePedestalRMSX;
      }
    }

#ifdef __DEBUG__
    float AdcSumBeforeAltro = 0, AdcSumAfterAltro = 0;
#endif

    for (int pad = 1; pad <= nPadsPerRow; pad++) {
      double gain = St_tpcPadGainT0BC::instance()->Gain(sector, row, pad);

      if (gain <= 0.0) continue;

      double ped = Cfg<TpcResponseSimulator>().AveragePedestal;
      static std::vector<short> ADCs(max_timebins_, 0);
      static std::vector<short> IDTs(max_timebins_, 0);
      std::fill(ADCs.begin(), ADCs.end(), 0);
      std::fill(IDTs.begin(), IDTs.end(), 0);
      int num_time_bins = 0;
      int index = max_timebins_ * ((row - 1) * max_pads_ + pad - 1);

      for (int bin = 0; bin < max_timebins_; bin++, index++) {
        int adc = pedRMS > 0 ? int(binned_charge[index].Sum / gain + gRandom->Gaus(ped, pedRMS)) - int(ped) :
                               int(binned_charge[index].Sum / gain);

        if (adc > 1023) adc = 1023;
        if (adc < 1) continue;

        num_time_bins++;
        ADCs[bin] = adc;
        IDTs[bin] = binned_charge[index].TrackId;
#ifdef __DEBUG__
        if (adc > 3 * pedRMS) AdcSumBeforeAltro += adc;
        if (Debug() > 11 && binned_charge[index].Sum > 0) {
          LOG_INFO << "digi R/P/T/I = " << row << " /\t" << pad << " /\t" << bin << " /\t" << index
                   << "\tSum/TrackId = " << binned_charge[index].Sum << " /\t" << binned_charge[index].TrackId << '\n';
        }
#endif
      }

      if (!num_time_bins) continue;

      if (Cfg<tpcAltroParams>(sector - 1).N >= 0) {
        Altro altro_sim(max_timebins_, ADCs.data());

        if (Cfg<tpcAltroParams>(sector - 1).N > 0) {
          //      ConfigAltro(ONBaselineCorrection1, ONTailcancellation, ONBaselineCorrection2, ONClipping, ONZerosuppression)
          altro_sim.ConfigAltro(                    0,                  1,                     0,          1,                 1);
          altro_sim.ConfigTailCancellationFilter(Cfg<tpcAltroParams>().Altro_K1,
                                              Cfg<tpcAltroParams>().Altro_K2,
                                              Cfg<tpcAltroParams>().Altro_K3,
                                              Cfg<tpcAltroParams>().Altro_L1,
                                              Cfg<tpcAltroParams>().Altro_L2,
                                              Cfg<tpcAltroParams>().Altro_L3);
        }
        else {
          altro_sim.ConfigAltro(0, 0, 0, 1, 1);
        }

        altro_sim.ConfigZerosuppression(Cfg<tpcAltroParams>().Altro_thr, Cfg<tpcAltroParams>().Altro_seq, 0, 0);
        altro_sim.RunEmulation();
        num_time_bins = 0;

        for (int i = 0; i < max_timebins_; i++) {
          if (ADCs[i] && !altro_sim.ADCkeep[i]) {ADCs[i] = 0;}

          if (ADCs[i]) {
            num_time_bins++;
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
      }
      else {
        if (Cfg<tpcAltroParams>(sector - 1).N < 0) num_time_bins = AsicThresholds(ADCs.data());
      }

      if (num_time_bins > 0) {
        digi_data.Add(sector, row, pad, ADCs.data(), IDTs.data(), max_timebins_);
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

  for (unsigned int tb = 0; tb < max_timebins_; tb++) {
    if (ADCs[tb] <= Cfg<asic_thresholds>().thresh_lo) {
      if (! t1) ADCs[tb] = 0;
      else {
        if (nSeqLo <= Cfg<asic_thresholds>().n_seq_lo ||
            nSeqHi <= Cfg<asic_thresholds>().n_seq_hi)
          for (unsigned int t = t1; t <= tb; t++) ADCs[t] = 0;
        else noTbleft += nSeqLo;
      }

      t1 = nSeqLo = nSeqHi = 0;
    }

    nSeqLo++;

    if (! t1) t1 = tb;

    if (ADCs[tb] > Cfg<asic_thresholds>().thresh_hi) {nSeqHi++;}
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
  // E.Mathieson (3.2b), V.Chernyatin said that it should be used this (Weber ) approximation 07/09/08
  double rc = s / (2 * M_PI) * std::exp(M_PI * h / s); LOG_INFO << "rc(Cylinder approx) = " << rc << " cm\n";
  //  double rc = 4*h/M_PI; LOG_INFO << "rc = " << rc << " cm\n";   // E.Mathieson (4.3), no valid for our case
  double C  = 1. / (2 * std::log(rc / ra)); LOG_INFO << "C = " << C << '\n';
  double E  = 2 * M_PI * C * Va / s; LOG_INFO << "E = " << E << " V/cm\n";
  // Gain variation: M = M0*(1 - k*cos(2*alpha))
  double k = 2 * B / 3.*std::pow((M_PI / E0 / s), 2) * std::pow(C * Va, 3); LOG_INFO << "k = " << k << '\n';
  // Induced charge variation
  t0 = ra * ra / (4 * mu * C * Va);
  LOG_INFO << "t0 = " << 1e9 * t0 << " ns\n";                                   // E.Mathieson (2.10)
  double Tav = t0 * h / s / (2 * M_PI * C);  LOG_INFO << "Tav = " << 1e9 * Tav << " ns\n";
  //  double t = 5*55e-9;             LOG_INFO << "t = " << 1e9*t << " ns\n";
  double t = 180e-9;             LOG_INFO << "t = " << 1e9 * t << " ns\n";
  double rp = std::sqrt(1. + t / t0); LOG_INFO << "r' = " << rp << '\n';
  // qc = rp*ra*sin(alpha)/(2*h) + C/2*log(1 + t/t0) = A*sin(alpha) + B
  double Aconstant = rp * ra / (2 * h);        LOG_INFO << "Aconstant = " << Aconstant << '\n';
  double Bconstant = C / 2 * std::log(1 + t / t0); LOG_INFO << "Bconstant = " << Bconstant << '\n';
  double Gains[2];

  for (int i = 0; i < 2; i++) {
    Gains[i] = Aconstant * std::sin(M_PI / 180 * alpha[i]) + Bconstant;
    LOG_INFO << "Gain = " << Gains[i] << " at alpha = " << alpha[i] << " degree\n";
  }

  double GainsAv = std::sqrt(Gains[0] * Gains[1]);
  double r = 0;

  for (int i = 0; i < 2; i++) {
    r = std::log(Gains[i] / GainsAv); LOG_INFO << "Relative gain " << r << " at alpha = " << alpha[i] << '\n';
  }

  return r;
}


double StTpcRSMaker::fei(double t, double t0, double T)
{
  static const double xmaxt = 708.39641853226408;
  static const double xmaxD  = xmaxt - std::log(xmaxt);
  double t01 = xmaxD, t11 = xmaxD;

  if (T > 0) {t11 = (t + t0) / T;}
  if (t11 > xmaxD) t11 = xmaxD;
  if (T > 0) {t01 = t0 / T;}
  if (t01 > xmaxD) t01  = xmaxD;

  return std::exp(-t11) * (ROOT::Math::expint(t11) - ROOT::Math::expint(t01));
}


double StTpcRSMaker::shapeEI(double* x, double* par)
{
  double t  = x[0];

  if (t <= 0) return 0;

  double t0    = par[0];
  double tau_I = par[3];
  double tau_C = par[5];

  if (std::abs((tau_I - tau_C) / (tau_I + tau_C)) < 1e-7) {
    return std::max(0., (t + t0) / tau_I * fei(t, t0, tau_I) + std::exp(-t / tau_I) - 1);
  }

  if (tau_C <= 0) return fei(t, t0, tau_I);
  if (tau_I <= 0) return 0;

  return tau_I / (tau_I - tau_C) * (fei(t, t0, tau_I) - fei(t, t0, tau_C));
}


double StTpcRSMaker::shapeEI3(double* x, double* par)
{
  double t  = x[0];
  double value = 0;

  if (t <= 0) return 0;

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
    value += A[i] * std::exp(a[i] * (t + t0)) * (ROOT::Math::expint(-a[i] * (t + t0)) - ROOT::Math::expint(-a[i] * t0));
  }

  return value;
}


double StTpcRSMaker::shapeEI_I(double* x, double* par)   //Integral of shape over time bin
{
  static double sqrt2 = std::sqrt(2.);
  double TimeBinWidth = par[4];
  double norm = par[5];
  double t1 = TimeBinWidth * (x[0] - 0.5);
  double t2 = t1 + TimeBinWidth;
  int io = (int) par[6];
  assert(io >= 0 && io <= 1);
  return sqrt2 * fgTimeShape0[io].Integral(t1, t2) / norm;
}


double StTpcRSMaker::shapeEI3_I(double* x, double* par)   //Integral of shape over time bin
{
  static double sqrt2 = std::sqrt(2.);
  double TimeBinWidth = par[4];
  double norm = par[5];
  double t1 = TimeBinWidth * (x[0] - 0.5);
  double t2 = t1 + TimeBinWidth;
  int io = (int) par[6];
  assert(io >= 0 && io <= 1);
  return sqrt2 * fgTimeShape3[io].Integral(t1, t2) / norm;
}


double StTpcRSMaker::polya(double* x, double* par)
{
  return tpcrs::GammaDist(x[0], par[0], par[1], par[2]);
}


double StTpcRSMaker::Ec(double* x, double* p)
{
  if (x[0] < p[0] / 2 || x[0] > 3.064 * p[0]) return 0;
  if (x[0] < p[0]) return 1;

  return std::pow(p[0] / x[0], 4);
}


void StTpcRSMaker::TrackSegment2Propagate(tpcrs::GeantHit& geant_hit, HitPoint_t &TrackSegmentHits, double& smin, double& smax)
{
  int volId = geant_hit.volume_id % 10000;
  int sector = volId / 100;

  TrackSegmentHits.xyzG = {geant_hit.x[0], geant_hit.x[1], geant_hit.x[2]};  PrPP(Make, TrackSegmentHits.xyzG);
  TrackSegmentHits.TrackId  = geant_hit.track_id;
  TrackSegmentHits.tpc_hitC = &geant_hit;
  TrackSegmentHits.sMin = TrackSegmentHits.s - TrackSegmentHits.tpc_hitC->ds;
  TrackSegmentHits.sMax = TrackSegmentHits.s;

  if (TrackSegmentHits.sMin < smin) smin = TrackSegmentHits.sMin;
  if (TrackSegmentHits.sMax > smax) smax = TrackSegmentHits.sMax;

  static StTpcLocalCoordinate coorLT;  // before do distortions
  static StTpcLocalSectorCoordinate coorS;
  static CoordTransform transform;
  // GlobalCoord -> LocalSectorCoord
  transform.global_to_local_sector(TrackSegmentHits.xyzG, coorS, sector, 0); PrPP(Make, coorS);
  int row = coorS.row;
  transform.global_to_local(TrackSegmentHits.xyzG, coorLT, sector, row); PrPP(Make, coorLT);

  // move up, calculate field at center of TPC
  static float BFieldG[3];
  StarMagField::Instance().BField(geant_hit.x, BFieldG);
  // distortion and misalignment
  // replace pxy => direction and try linear extrapolation
  Coords pxyzG{geant_hit.p[0], geant_hit.p[1], geant_hit.p[2]};
  StGlobalDirection     dirG{pxyzG.unit()};                                   PrPP(Make, dirG);
  StGlobalDirection     BG{BFieldG[0], BFieldG[1], BFieldG[2]};               PrPP(Make, BG);
  transform.global_to_local_sector_dir( dirG, TrackSegmentHits.dirLS, sector, row);  PrPP(Make, TrackSegmentHits.dirLS);
  transform.global_to_local_sector_dir(   BG, TrackSegmentHits.BLS,   sector, row);  PrPP(Make, TrackSegmentHits.BLS);

  // Distortions
  if (TESTBIT(options_, kDistortion) && StMagUtilities::Instance()) {
    float pos[3] = {(float ) coorLT.position.x, (float ) coorLT.position.y, (float ) coorLT.position.z};
    float posMoved[3];
    StMagUtilities::Instance()->DoDistortion(pos, posMoved, sector); // input pos[], returns posMoved[]
    coorLT.position = {posMoved[0], posMoved[1], posMoved[2]};        // after do distortions
    transform.local_to_global(coorLT, TrackSegmentHits.xyzG);                PrPP(Make, coorLT);
  }

  transform.local_to_local_sector(coorLT, TrackSegmentHits.coorLS); PrPP(Make, TrackSegmentHits.coorLS);

  double driftLength = TrackSegmentHits.coorLS.position.z + geant_hit.tof * StTpcDb::instance().DriftVelocity(sector); // ,row);

  if (driftLength > -1.0 && driftLength <= 0) {
    if ((row >  St_tpcPadConfigC::instance()->numberOfInnerRows(sector) && driftLength > - Cfg<tpcWirePlanes>().outerSectorAnodeWirePadSep) ||
        (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector) && driftLength > - Cfg<tpcWirePlanes>().innerSectorAnodeWirePadSep))
      driftLength = std::abs(driftLength);
  }

  TrackSegmentHits.coorLS.position.z = driftLength; PrPP(Make, TrackSegmentHits.coorLS);
  transform.local_sector_to_hardware(TrackSegmentHits.coorLS, TrackSegmentHits.Pad, false, false); // don't use T0, don't use Tau
  PrPP(Make, TrackSegmentHits.Pad);
}


std::vector<float> StTpcRSMaker::NumberOfElectronsInCluster(const TF1& heed, float dE, float& dEr)
{
  std::vector<float> rs;

  float dET = dE + dEr;
  dEr = dET;
  float EC;

  while ((EC = mHeed.GetRandom()) < dEr) {
    dEr -= EC;
    rs.push_back(1 - dEr / dET);
  }

  return rs;
}


double StTpcRSMaker::LoopOverElectronsInCluster(std::vector<float> rs, const HitPoint_t &TrackSegmentHits, std::vector<SignalSum_t>& binned_charge,
  int sector, int row,
  double xRange, Coords xyzC, double gain_local,
  double SigmaT, double SigmaL, double OmegaTau)
{
  double phiXY = 2 * M_PI * gRandom->Rndm();
  double rX = std::cos(phiXY);
  double rY = std::sin(phiXY);
  double total_signal_in_cluster = 0;
  int WireIndex = 0;

  int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;
  double sigmaJitterT = Cfg<TpcResponseSimulator>().SigmaJitterTI;
  double sigmaJitterX = Cfg<TpcResponseSimulator>().SigmaJitterXI;
  if (io) { // Outer
    sigmaJitterT = Cfg<TpcResponseSimulator>().SigmaJitterTO;
    sigmaJitterX = Cfg<TpcResponseSimulator>().SigmaJitterXO;
  }

  Coords unit = TrackSegmentHits.dirLS.position.unit();
  double L2L[9] = {unit.z,                  - unit.x*unit.z, unit.x,
                   unit.x,                  - unit.y*unit.z, unit.y,
                   0.0,       unit.x*unit.x + unit.y*unit.y, unit.z};

  for (int ie = 0; ie < rs.size(); ie++) {
    double gain_gas = mPolya[io].GetRandom();
    // transport to wire
    gRandom->Rannor(rX, rY);
    StTpcLocalSectorCoordinate xyzE{xyzC.x + rX * SigmaT,
                                    xyzC.y + rY * SigmaT,
                                    xyzC.z + gRandom->Gaus(0, SigmaL), sector, row};
    if (xRange > 0) {
      double xyzRangeL[3] = {rs[ie] * xRange * rX, rs[ie] * xRange * rY, 0.};
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
      for (int i=0; i<3; i++) xyzE.position[i] += xyzR[i];
    }

    double y = xyzE.position.y;
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
    else { // the first and last wires are fat ones
      WireIndex = tpcrs::irint((y - firstOuterSectorAnodeWire) / anodeWirePitch) + 1;

      if (WireIndex <= 1 || WireIndex >= numberOfOuterSectorAnodeWires) continue;

      yOnWire = firstOuterSectorAnodeWire + (WireIndex - 1) * anodeWirePitch;
      alphaVariation = OuterAlphaVariation[sector - 1];
    }

    double distanceToWire = y - yOnWire; // Calculated effective distance to wire affected by Lorentz shift
    xOnWire = xyzE.position.x;
    zOnWire = xyzE.position.z;
    // Grid plane (1 mm spacing) focusing effect + Lorentz angle in drift volume
    int iGridWire = int(std::abs(10.*distanceToWire));
    double dist2Grid = std::copysign(0.05 + 0.1 * iGridWire, distanceToWire); // [cm]
    // Ground plane (1 mm spacing) focusing effect
    int iGroundWire = int(std::abs(10.*dist2Grid));
    double distFocused = std::copysign(0.05 + 0.1 * iGroundWire, dist2Grid);
    // OmegaTau near wires taken from comparison with data
    double tanLorentz = OmegaTau / Cfg<TpcResponseSimulator>().OmegaTauScaleO;

    if (y < firstOuterSectorAnodeWire) tanLorentz = OmegaTau / Cfg<TpcResponseSimulator>().OmegaTauScaleI;

    xOnWire += distFocused * tanLorentz; // tanLorentz near wires taken from comparison with data
    zOnWire += std::abs(distFocused);

    if (! iGroundWire ) gain_gas *= std::exp( alphaVariation);
    else                gain_gas *= std::exp(-alphaVariation);

    double dY     = mChargeFraction[io][sector - 1].GetXmax();
    double yLmin  = yOnWire - dY;
    double yLmax  = yOnWire + dY;

    static CoordTransform transform;
    int    rowMin = transform.rowFromLocalY(yLmin, sector);
    int    rowMax = transform.rowFromLocalY(yLmax, sector);
    double yRmin  = transform.yFromRow(sector, rowMin) - St_tpcPadConfigC::instance()->PadLengthAtRow(sector, rowMin) / 2;
    double yRmax  = transform.yFromRow(sector, rowMax) + St_tpcPadConfigC::instance()->PadLengthAtRow(sector, rowMax) / 2;

    if (yRmin > yLmax || yRmax < yLmin) {
      continue;
    }

    GenerateSignal(TrackSegmentHits, sector, rowMin, rowMax, sigmaJitterT,
                   &mShaperResponses[io][sector - 1], binned_charge, total_signal_in_cluster, gain_local, gain_gas);
  }  // electrons in Cluster

  return total_signal_in_cluster;
}


void StTpcRSMaker::GenerateSignal(const HitPoint_t &TrackSegmentHits, int sector, int rowMin, int rowMax, double sigmaJitterT,
  TF1F* shaper, std::vector<SignalSum_t>& binned_charge, double& total_signal_in_cluster, double gain_local, double gain_gas)
{
  static CoordTransform transform;

  for (int row = rowMin; row <= rowMax; row++) {
    if (St_tpcPadConfigC::instance()->numberOfRows(sector) == 45) { // ! iTpx
      if ( !St_tpcRDOMasksC::instance()->isRowOn(sector, row) ) continue;
      if ( !St_tpcAnodeHVavgC::instance()->livePadrow(sector, row) )  continue;
    }

    int io = (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ? 0 : 1;
    StTpcLocalSectorCoordinate xyzW{xOnWire, yOnWire, zOnWire, sector, row};
    static StTpcPadCoordinate Pad;
    transform.local_sector_to_hardware(xyzW, Pad, false, false); // don't use T0, don't use Tau
    float bin = Pad.timeBucket;//L  - 1; // K
    int binT = tpcrs::irint(bin); //L bin;//K tpcrs::irint(bin);// J bin; // I tpcrs::irint(bin);

    if (binT < 0 || binT >= max_timebins_) continue;

    double dT = bin - binT + Cfg<TpcResponseSimulator>().T0offset;
    dT += (row <= St_tpcPadConfigC::instance()->numberOfInnerRows(sector)) ?
          Cfg<TpcResponseSimulator>().T0offsetI :
          Cfg<TpcResponseSimulator>().T0offsetO;

    if (sigmaJitterT) dT += gRandom->Gaus(0, sigmaJitterT); // #1

    double dely = transform.yFromRow(sector, row) - yOnWire;
    double localYDirectionCoupling = mChargeFraction[io][sector - 1].GetSaveL(&dely);

    if (localYDirectionCoupling < min_signal_) continue;

    float padX = Pad.pad;
    int CentralPad = tpcrs::irint(padX);
    int PadsAtRow = St_tpcPadConfigC::instance()->numberOfPadsAtRow(sector, row);

    if (CentralPad < 1 || CentralPad > PadsAtRow) continue;

    int DeltaPad = tpcrs::irint(mPadResponseFunction[io][sector - 1].GetXmax()) + 1;
    int padMin   = std::max(CentralPad - DeltaPad, 1);
    int padMax   = std::min(CentralPad + DeltaPad, PadsAtRow);
    int Npads    = std::min(padMax - padMin + 1, static_cast<int>(kPadMax));
    double xPadMin = padMin - padX;
    static double XDirectionCouplings[kPadMax];
    static double TimeCouplings[kTimeBacketMax];
    mPadResponseFunction[io][sector - 1].GetSaveL(Npads, xPadMin, XDirectionCouplings);

    for (int pad = padMin; pad <= padMax; pad++) {
      double gain = gain_gas * gain_local;
      double dt = dT;

      if (! TESTBIT(options_, kGAINOAtALL)) {
        gain *= St_tpcPadGainT0BC::instance()->Gain(sector, row, pad);

        if (gain <= 0.0) continue;

        dt -= St_tpcPadGainT0BC::instance()->T0(sector, row, pad);
      }

      double localXDirectionCoupling = gain * XDirectionCouplings[pad - padMin];

      if (localXDirectionCoupling < min_signal_) continue;

      double XYcoupling = localYDirectionCoupling * localXDirectionCoupling;

      if (XYcoupling < min_signal_)  continue;

      int bin_low  = std::max(0, binT + tpcrs::irint(dt + shaper->GetXmin() - 0.5));
      int bin_high = std::min(max_timebins_ - 1, binT + tpcrs::irint(dt + shaper->GetXmax() + 0.5));
      int index = max_timebins_ * ((row - 1) * max_pads_ + pad - 1) + bin_low;
      int Ntbks = std::min(bin_high - bin_low + 1, static_cast<int>(kTimeBacketMax));
      double tt = -dt + (bin_low - binT);
      shaper->GetSaveL(Ntbks, tt, TimeCouplings);

      for (int itbin = bin_low; itbin <= bin_high; itbin++, index++) {
        double signal = XYcoupling * TimeCouplings[itbin - bin_low];

        if (signal < min_signal_)  continue;

        total_signal_in_cluster += signal;
        binned_charge[index].Sum += signal;

        rowsdE[row - 1]  += signal;

        if ( TrackSegmentHits.TrackId ) {
          if ( !binned_charge[index].TrackId )
            binned_charge[index].TrackId = TrackSegmentHits.TrackId;
          else { // switch TrackId, works only for 2 tracks, more tracks ?
            if ( binned_charge[index].TrackId != TrackSegmentHits.TrackId && binned_charge[index].Sum < 2 * signal)
              binned_charge[index].TrackId = TrackSegmentHits.TrackId;
          }
        }

#ifdef __DEBUG__
        if (Debug() > 13 && (binned_charge[index].Sum > 0 || ! std::isfinite(binned_charge[index].Sum)) ) {
          LOG_INFO << "simu row = " << TrackSegmentHits.tpc_hitC->volume_id % 100 << "\tR/P/T/I = " << row << " /\t" << pad << " /\t" << itbin << " /\t" << index
                   << "\tSum/TrackId = " << binned_charge[index].Sum << " /\t" << binned_charge[index].TrackId
                   << "\tsignal = " << signal
                   << "\trow Min/Max = " << rowMin << "/" << rowMax
                   << '\n';

          if (! std::isfinite(binned_charge[index].Sum)) {
            LOG_INFO << "Not Finite\n";
          }
        }
#endif
      } // time
    } // pad limits
  } // row limits
}


double StTpcRSMaker::dEdxCorrection(const HitPoint_t &path_segment)
{
  double dEdxCor = 1;
  double dStep =  std::abs(path_segment.tpc_hitC->ds);
  dEdxY2_t CdEdx;
  memset (&CdEdx, 0, sizeof(dEdxY2_t));
  CdEdx.DeltaZ = 5.2;
  CdEdx.QRatio = -2;
  CdEdx.QRatioA = -2.;
  CdEdx.QSumA = 0;
  CdEdx.sector = path_segment.Pad.sector;
  CdEdx.row    = path_segment.Pad.row;
  CdEdx.pad    = tpcrs::irint(path_segment.Pad.pad);
  CdEdx.edge   = CdEdx.pad;

  if (CdEdx.edge > 0.5 * St_tpcPadConfigC::instance()->numberOfPadsAtRow(CdEdx.sector, CdEdx.row))
    CdEdx.edge += 1 - St_tpcPadConfigC::instance()->numberOfPadsAtRow(CdEdx.sector, CdEdx.row);

  CdEdx.F.dE     = 1;
  CdEdx.F.dx     = dStep;
  CdEdx.xyz[0] = path_segment.coorLS.position.x;
  CdEdx.xyz[1] = path_segment.coorLS.position.y;
  CdEdx.xyz[2] = path_segment.coorLS.position.z;
  double probablePad = St_tpcPadConfigC::instance()->numberOfPadsAtRow(CdEdx.sector, CdEdx.row) / 2;
  double pitch = (CdEdx.row <= St_tpcPadConfigC::instance()->numberOfInnerRows(CdEdx.sector)) ?
                   St_tpcPadConfigC::instance()->innerSectorPadPitch(CdEdx.sector) :
                   St_tpcPadConfigC::instance()->outerSectorPadPitch(CdEdx.sector);
  double PhiMax = std::atan2(probablePad * pitch, St_tpcPadConfigC::instance()->radialDistanceAtRow(CdEdx.sector, CdEdx.row));
  CdEdx.PhiR   = std::atan2(CdEdx.xyz[0], CdEdx.xyz[1]) / PhiMax;
  CdEdx.xyzD[0] = path_segment.dirLS.position.x;
  CdEdx.xyzD[1] = path_segment.dirLS.position.y;
  CdEdx.xyzD[2] = path_segment.dirLS.position.z;
  CdEdx.ZdriftDistance = CdEdx.xyzD[2];
  CdEdx.zG      = CdEdx.xyz[2];

  if (St_trigDetSumsC::instance())	CdEdx.Zdc     = St_trigDetSumsC::instance()->zdcX();

  CdEdx.ZdriftDistance = path_segment.coorLS.position.z; // drift length
  St_tpcGasC* tpc_gas = m_TpcdEdxCorrection.TpcGas();

  if (tpc_gas)
    CdEdx.ZdriftDistanceO2 = CdEdx.ZdriftDistance * tpc_gas->Struct()->ppmOxygenIn;

  if (! m_TpcdEdxCorrection.dEdxCorrection(CdEdx)) {
    dEdxCor = CdEdx.F.dE;
  }

  return dEdxCor;
}

#undef PrPP
