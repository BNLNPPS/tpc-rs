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
#include "mag_field_utils.h"
#include "math_funcs.h"
#include "struct_containers.h"
#include "track_helix.h"


#define __DEBUG__
#if defined(__DEBUG__)
#define PrPP(A,B) if (Debug()%10 > 2) {LOG_INFO << "Simulator::" << (#A) << "\t" << (#B) << " = \t" << (B) << '\n';}
#else
#define PrPP(A,B)
#endif

//                                    Inner        Outer
static       double t0IO[2]   = {1.20868e-9, 1.43615e-9}; // recalculated in InducedCharge

TF1F Simulator::fgTimeShape3[2] = {
  TF1F("TimeShape3Inner;Time [s];Signal", Simulator::shapeEI3, 0, 1, 7),
  TF1F("TimeShape3Outer;Time [s];Signal", Simulator::shapeEI3, 0, 1, 7)
};

TF1F Simulator::fgTimeShape0[2] = {
  TF1F("TimeShape0Inner;Time [s];Signal", Simulator::shapeEI, 0, 1, 7),
  TF1F("TimeShape0Outer;Time [s];Signal", Simulator::shapeEI, 0, 1, 7)
};

using dEdxCorr = StTpcdEdxCorrection::Corrections;

Simulator::Simulator(const tpcrs::Configurator& cfg):
  cfg_(cfg),
  transform_(cfg),
  mag_field_utils_(cfg, transform_),
  num_sectors_(cfg_.S<tpcDimensions>().numberOfSectors),
  max_rows_(72),
  max_pads_(182),
  max_timebins_(512),
  options_(0),
  mdNdx(nullptr),
  mdNdxL10(nullptr),
  mdNdEL10(nullptr),
  mShaperResponses{
    std::vector<TF1F>(num_sectors_, TF1F("ShaperFuncInner;Time [bin];Signal", Simulator::shapeEI_I, 0, 1, 7)),
    std::vector<TF1F>(num_sectors_, TF1F("ShaperFuncOuter;Time [bin];Signal", Simulator::shapeEI_I, 0, 1, 7))
  },
  mChargeFraction{
    std::vector<TF1F>(num_sectors_, TF1F("ChargeFractionInner;Distance [cm];Signal", Simulator::PadResponseFunc, 0, 1, 6)),
    std::vector<TF1F>(num_sectors_, TF1F("ChargeFractionOuter;Distance [cm];Signal", Simulator::PadResponseFunc, 0, 1, 6))
  },
  mPadResponseFunction{
    std::vector<TF1F>(num_sectors_, TF1F("PadResponseFunctionInner;Distance [pads];Signal", Simulator::PadResponseFunc, 0, 1, 6)),
    std::vector<TF1F>(num_sectors_, TF1F("PadResponseFunctionOuter;Distance [pads];Signal", Simulator::PadResponseFunc, 0, 1, 6))
  },
  mPolya{
    TF1F("PolyaInner;x = G/G_0;signal", polya, 0, 10, 3),
    TF1F("PolyaOuter;x = G/G_0;signal", polya, 0, 10, 3)
  },
  mHeed("Ec", Simulator::Ec, 0, 3.064 * cfg_.S<TpcResponseSimulator>().W, 1),
  dEdx_correction_(cfg, dEdxCorr::kAll & ~dEdxCorr::kAdcCorrection & ~dEdxCorr::kAdcCorrectionMDF & ~dEdxCorr::kdXCorrection, Debug())
{
  //  SETBIT(options_,kHEED);
  SETBIT(options_, kBICHSEL); // Default is Bichsel
  SETBIT(options_, kdEdxCorr);
  //SETBIT(options_, kDistortion);

  if (TESTBIT(options_, kBICHSEL)) {
    LOG_INFO << "Simulator:: use H.Bichsel model for dE/dx simulation\n";

    TFile inner(cfg_.Locate("dNdE_Bichsel.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(cfg_.Locate("dNdx_Bichsel.root").c_str());
    mdNdx = (TH1D*) outer.Get("dNdx"); assert(mdNdx);
    mdNdx->SetDirectory(0);
  }
  else if (TESTBIT(options_, kHEED)) {
    LOG_INFO << "Simulator:: use Heed model for dE/dx simulation\n";

    TFile inner(cfg_.Locate("dNdx_Heed.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(cfg_.Locate("dNdx_Heed.root").c_str());
    mdNdxL10 = (TH1D*) outer.Get("dNdxL10"); assert(mdNdxL10);
    mdNdxL10->SetDirectory(0);
  }
  else {LOG_INFO << "Simulator:: use GEANT321 model for dE/dx simulation\n";}

  double TimeBinWidth = 1. / cfg_.S<starClockOnl>().frequency;
  numberOfInnerSectorAnodeWires  = cfg_.S<tpcWirePlanes>().numInnerSectorAnodeWires;
  firstInnerSectorAnodeWire      = cfg_.S<tpcWirePlanes>().firstInnerSectorAnodeWire;
  lastInnerSectorAnodeWire       = cfg_.S<tpcWirePlanes>().lastInnerSectorAnodeWire;
  numberOfOuterSectorAnodeWires  = cfg_.S<tpcWirePlanes>().numOuterSectorAnodeWires;
  firstOuterSectorAnodeWire      = cfg_.S<tpcWirePlanes>().firstOuterSectorAnodeWire;
  lastOuterSectorAnodeWire       = cfg_.S<tpcWirePlanes>().lastOuterSectorAnodeWire;
  anodeWirePitch                 = cfg_.S<tpcWirePlanes>().anodeWirePitch;
  anodeWireRadius                = cfg_.S<tpcWirePlanes>().anodeWireRadius;

  // Shapers
  double timeBinMin = -0.5;
  double timeBinMax = 44.5;
  double CathodeAnodeGap[2] = {0.2, 0.4};

  for (int sector = 1; sector <= num_sectors_; sector++) {
    innerSectorAnodeVoltage[sector - 1] = outerSectorAnodeVoltage[sector - 1] = 0;
    int nAliveInner = 0;
    int nAliveOuter = 0;

    for (int row = 1; row <= cfg_.C<St_tpcPadConfigC>().numberOfRows(sector); row++) {
      if (cfg_.C<St_tpcPadConfigC>().IsRowInner(sector, row)) {
        nAliveInner++;
        innerSectorAnodeVoltage[sector - 1] += cfg_.C<St_tpcAnodeHVavgC>().voltagePadrow(sector, row);
      }
      else {
        nAliveOuter++;
        outerSectorAnodeVoltage[sector - 1] += cfg_.C<St_tpcAnodeHVavgC>().voltagePadrow(sector, row);
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
      {"tauF",  cfg_.S<TpcResponseSimulator>().tauF},
      {"tauP",  cfg_.S<TpcResponseSimulator>().tauP},
      {"tauI",  cfg_.S<TpcResponseSimulator>().tauIntegration},
      {"width", TimeBinWidth},
      {"tauC",  0},
      {"io",    io}
    };

    FuncParams_t params0{
      {"t0",    t0IO[io]},
      {"tauF",  0},
      {"tauP",  0},
      {"tauI",  io ? cfg_.S<TpcResponseSimulator>().tauXO : cfg_.S<TpcResponseSimulator>().tauXI},
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
    params0[5].second = io ? cfg_.S<TpcResponseSimulator>().tauCO : cfg_.S<TpcResponseSimulator>().tauCI;
    for (int i = 0; i != params0.size(); ++i) {
      fgTimeShape0[io].SetParName(i, params0[i].first.c_str());
      fgTimeShape0[io].SetParameter(i, params0[i].second);
    }
    fgTimeShape0[io].SetRange(timeBinMin * TimeBinWidth, timeBinMax * TimeBinWidth);
    params0[5].second = fgTimeShape0[io].Integral(0, timeBinMax * TimeBinWidth);

    for (int sector = 1; sector <= num_sectors_; sector++) {
      //                            w       h        s       a       l  i
      //  double paramsI[6] = {0.2850, 0.2000,  0.4000, 0.0010, 1.1500, 0};
      //  double paramsO[6] = {0.6200, 0.4000,  0.4000, 0.0010, 1.1500, 0};
      double params[6]{
        io == kInner ? cfg_.C<St_tpcPadConfigC>().innerSectorPadWidth(sector) :  // w = width of pad
                       cfg_.C<St_tpcPadConfigC>().outerSectorPadWidth(sector),
        io == kInner ? cfg_.S<tpcWirePlanes>().innerSectorAnodeWirePadSep :            // h = Anode-Cathode gap
                       cfg_.S<tpcWirePlanes>().outerSectorAnodeWirePadSep,
        cfg_.S<tpcWirePlanes>().anodeWirePitch,                                        // s = wire spacing
        io == kInner ? cfg_.S<TpcResponseSimulator>().K3IP :
                       cfg_.S<TpcResponseSimulator>().K3OP,
        0,
        io == kInner ? cfg_.C<St_tpcPadConfigC>().innerSectorPadPitch(sector) :
                       cfg_.C<St_tpcPadConfigC>().outerSectorPadPitch(sector)
      };

      mPadResponseFunction[io][sector - 1].SetParameters(params);
      mPadResponseFunction[io][sector - 1].SetParNames("PadWidth", "Anode-Cathode gap", "wire spacing", "K3OP", "CrossTalk", "PadPitch");
      mPadResponseFunction[io][sector - 1].SetRange(-2.5, 2.5); // Cut tails
      mPadResponseFunction[io][sector - 1].Save(-4.5, 4.5, 0, 0, 0, 0);

      params[0] = io == kInner ? cfg_.C<St_tpcPadConfigC>().innerSectorPadLength(sector) :
                                 cfg_.C<St_tpcPadConfigC>().outerSectorPadLength(sector);
      params[3] = io == kInner ? cfg_.S<TpcResponseSimulator>().K3IR :
                                 cfg_.S<TpcResponseSimulator>().K3OR;
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
      if (cfg_.S<tpcAltroParams>(sector - 1).N < 0) { // old TPC
        InitShaperFuncs(io, sector, mShaperResponses, Simulator::shapeEI3_I, params3, timeBinMin, timeBinMax);
      } else {//Altro
        InitShaperFuncs(io, sector, mShaperResponses, Simulator::shapeEI_I,  params0, timeBinMin, timeBinMax);
      }
    }
  }

  if (Debug()) Print();

  //  mPolya = new TF1F("Polya;x = G/G_0;signal","sqrt(x)/exp(1.5*x)",0,10); // original Polya
  //  mPolya = new TF1F("Polya;x = G/G_0;signal","pow(x,0.38)*exp(-1.38*x)",0,10); //  Valeri Cherniatin
  //   mPoly = new TH1D("Poly","polyaAvalanche",100,0,10);
  //if (gamma <= 0) gamma = 1.38;
  double gamma_inn = cfg_.S<TpcResponseSimulator>().PolyaInner;
  double gamma_out = cfg_.S<TpcResponseSimulator>().PolyaOuter;
  mPolya[kInner].SetParameters(gamma_inn, 0., 1. / gamma_inn);
  mPolya[kOuter].SetParameters(gamma_out, 0., 1. / gamma_out);

  // HEED function to generate Ec, default w = 26.2
  mHeed.SetParameter(0, cfg_.S<TpcResponseSimulator>().W);
}


Simulator::~Simulator()
{
  delete mdNdx;
  delete mdNdxL10;
  delete mdNdEL10;
}


void Simulator::InitShaperFuncs(int io, int sector, std::array<std::vector<TF1F>, 2>& funcs,
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


template<>
void Simulator::Simulate(GeneratedHitIt first_hit, GeneratedHitIt last_hit, DigiInserter digi_data)
{
  static int nCalls = 0;
  gRandom->SetSeed(2345 + nCalls++);

  double vminI = cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min;
  double vminO = cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min;

  // TODO: Confirm proper handling of empty input containers

  cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min = -500;
  cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min = -500;

  if (Debug()) {
    LOG_INFO << "Reset min for gain Correction to I/O\t"
             << cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min
             << "\t"
             << cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min
             << " (V)\n";
  }

  std::vector< std::vector<TrackSegment> > segments_by_sector(num_sectors_);
  std::vector<TrackSegment> segments_in_sector;

  auto first_hit_on_track = first_hit;
  int curr_direction = 0; // 0 - increase no of row, 1 - decrease no of. row.

  for (auto curr_hit = first_hit; curr_hit != last_hit; ++curr_hit)
  {
    auto next_hit = next(curr_hit);    

    int  next_direction  =  next_hit->volume_id % 100 - curr_hit->volume_id % 100 >= 0 ? 0 : 1;
    bool sector_boundary = (next_hit->volume_id % 10000) / 100 !=
                           (curr_hit->volume_id % 10000) / 100;
    bool track_boundary  = (next_hit->track_id != curr_hit->track_id || curr_hit->track_id == 0);

    bool start_new_track = track_boundary || curr_direction != next_direction || sector_boundary;

    if (start_new_track || next_hit == last_hit) {
      CreateTrackSegments(first_hit_on_track, next_hit, segments_in_sector);
      first_hit_on_track = next_hit;

      if ( (sector_boundary || next_hit == last_hit ) && segments_in_sector.size() != 0) {
        segments_by_sector[(curr_hit->volume_id % 10000) / 100 - 1] = segments_in_sector;
        segments_in_sector.clear();
      }
    }

    curr_direction = next_direction;
  }

  int sector = 1;
  for (auto segments_in_sector : segments_by_sector) {
    int nHitsInTheSector = 0;
    ChargeContainer binned_charge(max_rows_ * max_pads_ * max_timebins_);

    for (TrackSegment& segment : segments_in_sector) {
      segment.tpc_hitC->digi.adc = 0;
      int row = segment.coorLS.row;

      double gain_base = CalcBaseGain(segment.Pad.sector, segment.Pad.row);

      // dE/dx correction
      double dEdxCor = dEdxCorrection(segment);
      dEdxCor *= GatingGridTransparency(segment.Pad.timeBucket);
      if (dEdxCor < cfg_.S<ResponseSimulator>().min_signal) continue;

      double gain_local = CalcLocalGain(segment.Pad.sector, segment.Pad.row, gain_base, dEdxCor);

      // Initialize propagation
      // Magnetic field BField must be in kilogauss
      // kilogauss = 1e-1*tesla = 1e-1*(volt*second/meter2) = 1e-1*(1e-6*1e-3*1/1e4) = 1e-14
      TrackHelix track(segment.dirLS.position,
                       segment.coorLS.position,
                       segment.BLS.position.z * 1e-14 * segment.charge, 1);
      // Propagate track to middle of the pad row plane by the nominal center point and the normal
      // in this sector coordinate system
      double sR = track.pathLength({0, transform_.yFromRow(segment.Pad.sector, segment.Pad.row), 0}, {0, 1, 0});

      // Update hit position based on the new track crossing the middle of pad row
      if (sR < 1e10) {
        PrPP(Make, sR);
        PrPP(Make, segment.coorLS);
        segment.coorLS.position = {track.at(sR).x, track.at(sR).y, track.at(sR).z};
        PrPP(Make, segment.coorLS);
        PrPP(Make, segment.Pad);
        transform_.local_sector_to_hardware(segment.coorLS, segment.Pad, false, false); // don't use T0, don't use Tau
        PrPP(Make, segment.Pad);
      }

      int nP = 0;
      double dESum = 0;
      double dSSum = 0;

      SignalFromSegment(segment, track, gain_local, binned_charge, nP, dESum, dSSum);

      static const double eV = 1e-9; // electronvolt in GeV
      segment.tpc_hitC->digi.de = dESum * eV;
      segment.tpc_hitC->digi.ds = dSSum;
      segment.tpc_hitC->digi.np = nP;
      segment.tpc_hitC->digi.pad = segment.Pad.pad;
      segment.tpc_hitC->digi.timebin = segment.Pad.timeBucket;

      nHitsInTheSector++;
    } // end do loop over segments for a given particle

    if (nHitsInTheSector) {
      DigitizeSector(sector, digi_data, binned_charge);

      if (Debug()) LOG_INFO << "Simulator: Done with sector\t" << sector << " total no. of hit = " << nHitsInTheSector << '\n';
    }
    sector++;
  } // sector

  cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min = vminI;
  cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min = vminO;

  if (Debug()) {
    LOG_INFO << "Reset min for gain Correction to I/O\t"
             << cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min
             << "\t"
             << cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min
             << " (V)\n";
  }
}


void Simulator::CreateTrackSegments(GeneratedHitIt first_hit, GeneratedHitIt last_hit, std::vector<TrackSegment>& segments)
{
  for (auto ihit = first_hit; ihit != last_hit; ++ihit)
  {
    TrackSegment curr_segment = CreateTrackSegment(*ihit);

    if (curr_segment.charge == 0) continue;
    if (curr_segment.Pad.timeBucket < 0 || curr_segment.Pad.timeBucket > max_timebins_) continue;

    segments.push_back(curr_segment);
  }
}


double Simulator::GetNoPrimaryClusters(double betaGamma, int charge)
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


double Simulator::PadResponseFunc(double* x, double* par)
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


double Simulator::Gatti(double* x, double* par)
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


void  Simulator::Print(Option_t* /* option */) const
{
  PrPP(Print, num_sectors_);
  PrPP(Print, cfg_.C<St_tpcPadConfigC>().numberOfRows(1));
  PrPP(Print, cfg_.C<St_tpcPadConfigC>().numberOfRows(20));
  PrPP(Print, cfg_.C<St_tpcPadConfigC>().numberOfInnerRows(20));
  PrPP(Print, max_pads_);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().W);// = 26.2);//*eV
  PrPP(Print, cfg_.S<TpcResponseSimulator>().Cluster);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().longitudinalDiffusion);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().transverseDiffusion);
  //  PrPP(Print, Gain);
  PrPP(Print, max_timebins_);
  PrPP(Print, numberOfInnerSectorAnodeWires);
  PrPP(Print, firstInnerSectorAnodeWire);
  PrPP(Print, lastInnerSectorAnodeWire);
  PrPP(Print, numberOfOuterSectorAnodeWires);
  PrPP(Print, firstOuterSectorAnodeWire);
  PrPP(Print, lastOuterSectorAnodeWire);
  PrPP(Print, anodeWirePitch);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().OmegaTau); // tan of Lorentz angle
  PrPP(Print, cfg_.S<TpcResponseSimulator>().NoElPerAdcI);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().NoElPerAdcO);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().NoElPerAdcX);
  PrPP(Print, anodeWireRadius);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().AveragePedestal);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().AveragePedestalRMS);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().FanoFactor);

  for (int sector = 1; sector <= num_sectors_; sector++) {
    PrPP(Print, innerSectorAnodeVoltage[sector-1]);
    PrPP(Print, outerSectorAnodeVoltage[sector-1]);
  }

  PrPP(Print, cfg_.S<TpcResponseSimulator>().K3IP);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().K3IR);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().K3OP);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().K3OR);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().SigmaJitterTI);
  PrPP(Print, cfg_.S<TpcResponseSimulator>().SigmaJitterTO);
}


void Simulator::DigitizeSector(int sector, DigiInserter digi_data, const ChargeContainer& binned_charge)
{
  for (int row = 1;  row <= cfg_.C<St_tpcPadConfigC>().numberOfRows(sector); row++) {
    int nPadsPerRow = cfg_.C<St_tpcPadConfigC>().padsPerRow(sector, row);
    double pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMS;

    if (cfg_.S<tpcAltroParams>(sector - 1).N > 0) {
      if (! (cfg_.C<St_tpcPadConfigC>().iTPC(sector) && cfg_.C<St_tpcPadConfigC>().IsRowInner(sector, row))) {
        pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX;
      }
    }

    for (int pad = 1; pad <= nPadsPerRow; pad++) {
      double gain = cfg_.C<St_tpcPadGainT0BC>().Gain(sector, row, pad);

      if (gain <= 0.0) continue;

      double ped = cfg_.S<TpcResponseSimulator>().AveragePedestal;
      static std::vector<short> ADCs(max_timebins_, 0);
      static std::vector<short> IDTs(max_timebins_, 0);
      std::fill(ADCs.begin(), ADCs.end(), 0);
      std::fill(IDTs.begin(), IDTs.end(), 0);
      int index = max_timebins_ * ((row - 1) * max_pads_ + pad - 1);

      for (int bin = 0; bin < max_timebins_; bin++, index++) {
        int adc = pedRMS > 0 ? int(binned_charge[index].Sum / gain + gRandom->Gaus(ped, pedRMS)) - int(ped) :
                               int(binned_charge[index].Sum / gain);

        if (adc > 1023) adc = 1023;
        if (adc < 1) continue;

        ADCs[bin] = adc;
        IDTs[bin] = binned_charge[index].TrackId;
      }

      int flag = cfg_.S<tpcAltroParams>(sector - 1).N;
      flag >= 0 ? SimulateAltro(ADCs.begin(), ADCs.end(), flag>0) : SimulateAsic(ADCs);

      AddDigiData(sector, row, pad, ADCs.data(), IDTs.data(), max_timebins_, digi_data);
    } // pads
  } // row
}


void Simulator::AddDigiData(unsigned int sector, unsigned int row, unsigned int pad, short* ADCs, short* IDTs, int n_timebins, DigiInserter digi_data)
{
  bool in_cluster = false;

  for (unsigned int tb = 0; tb < n_timebins; ++tb)
  {
    if (!ADCs[tb])
      in_cluster = false;

    if (ADCs[tb] && !in_cluster)
      in_cluster = true;

    if (in_cluster)
      *digi_data = tpcrs::DigiHit{sector, row, pad, tb, ADCs[tb], IDTs[tb]};
  }
}


void Simulator::SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail)
{
  Altro altro_sim(last - first, &*first);

  if (cancel_tail) {
    //        ConfigAltro(ONBaselineCorrection1, ONTailcancellation, ONBaselineCorrection2, ONClipping, ONZerosuppression)
    altro_sim.ConfigAltro(                    0,                  1,                     0,          1,                 1);
    altro_sim.ConfigTailCancellationFilter(cfg_.S<tpcAltroParams>().Altro_K1,
                                        cfg_.S<tpcAltroParams>().Altro_K2,
                                        cfg_.S<tpcAltroParams>().Altro_K3,
                                        cfg_.S<tpcAltroParams>().Altro_L1,
                                        cfg_.S<tpcAltroParams>().Altro_L2,
                                        cfg_.S<tpcAltroParams>().Altro_L3);
  }
  else {
    altro_sim.ConfigAltro(0, 0, 0, 1, 1);
  }

  altro_sim.ConfigZerosuppression(cfg_.S<tpcAltroParams>().Altro_thr, cfg_.S<tpcAltroParams>().Altro_seq, 0, 0);
  altro_sim.RunEmulation();

  int i = 0;
  for (auto it = first; it != last; ++it, ++i) {
    if (*it && !altro_sim.ADCkeep[i]) { *it = 0;}
  }
}


void Simulator::SimulateAsic(std::vector<short>& ADC)
{
  int t1 = 0;
  int nSeqLo = 0;
  int nSeqHi = 0;

  for (unsigned int tb = 0; tb != ADC.size(); ++tb) {
    if (ADC[tb] <= cfg_.S<asic_thresholds>().thresh_lo) {
      if (! t1) ADC[tb] = 0;
      else {
        if (nSeqLo <= cfg_.S<asic_thresholds>().n_seq_lo ||
            nSeqHi <= cfg_.S<asic_thresholds>().n_seq_hi)
        {
          for (unsigned int t = t1; t <= tb; t++) ADC[t] = 0;
        }
      }

      t1 = nSeqLo = nSeqHi = 0;
    }

    nSeqLo++;

    if (!t1) t1 = tb;
    if (ADC[tb] > cfg_.S<asic_thresholds>().thresh_hi) {nSeqHi++;}
  }
}


double Simulator::InducedCharge(double s, double h, double ra, double Va, double &t0)
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


void Simulator::ParticleProperties(int particle_id, int& charge, double& mass)
{
  const int LASERINO = 170;
  const int CHASRINO = 171;

  charge = 0;
  mass = 0;

  StParticleDefinition* particle = StParticleTable::instance()->findParticleByGeantId(particle_id);

  if (particle) {
    mass = particle->mass();
    charge = particle->charge();
  }

  if (particle_id == LASERINO || particle_id == CHASRINO) {
    charge = 0;
  }
  else {
    if (particle_id == 1) {// gamma => electron
      particle_id = 3;
      charge = -1;
    }
  }
  // special treatment for electron/positron
  if (particle_id == 2) charge =  101;
  if (particle_id == 3) charge = -101;
};


double Simulator::fei(double t, double t0, double t1)
{
  static const double xmaxt = 708.39641853226408;
  static const double xmaxD  = xmaxt - std::log(xmaxt);
  double t01 = xmaxD, t11 = xmaxD;

  if (t1 > 0) {t11 = (t + t0) / t1;}
  if (t11 > xmaxD) t11 = xmaxD;
  if (t1 > 0) {t01 = t0 / t1;}
  if (t01 > xmaxD) t01  = xmaxD;

  return std::exp(-t11) * (ROOT::Math::expint(t11) - ROOT::Math::expint(t01));
}


double Simulator::shapeEI(double* x, double* par)
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


double Simulator::shapeEI3(double* x, double* par)
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


double Simulator::shapeEI_I(double* x, double* par)   //Integral of shape over time bin
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


double Simulator::shapeEI3_I(double* x, double* par)   //Integral of shape over time bin
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


double Simulator::polya(double* x, double* par)
{
  return tpcrs::GammaDist(x[0], par[0], par[1], par[2]);
}


double Simulator::Ec(double* x, double* p)
{
  if (x[0] < p[0] / 2 || x[0] > 3.064 * p[0]) return 0;
  if (x[0] < p[0]) return 1;

  return std::pow(p[0] / x[0], 4);
}


Simulator::TrackSegment Simulator::CreateTrackSegment(tpcrs::SimulatedHit& hit)
{
  int volId = hit.volume_id % 10000;
  int sector = volId / 100;

  TrackSegment segment;
  StGlobalCoordinate xyzG{hit.x[0], hit.x[1], hit.x[2]};  PrPP(Make, xyzG);
  segment.tpc_hitC = &hit;
  ParticleProperties(hit.particle_id, segment.charge, segment.mass);

  StTpcLocalSectorCoordinate coorS;
  // GlobalCoord -> LocalSectorCoord. This transformation can result in a row
  // that is not the same as (volId % 100)
  transform_.global_to_local_sector(xyzG, coorS, sector, 0); PrPP(Make, coorS);
  StTpcLocalCoordinate coorLT;  // before distortions
  transform_.global_to_local(xyzG, coorLT, sector, coorS.row); PrPP(Make, coorLT);

  // move up, calculate field at center of TPC
  static float BFieldG[3];
  mag_field_utils_.BFieldTpc(hit.x, BFieldG);
  // distortion and misalignment
  // replace pxy => direction and try linear extrapolation
  Coords pxyzG{hit.p[0], hit.p[1], hit.p[2]};
  StGlobalDirection dirG{pxyzG.unit()};                                      PrPP(Make, dirG);
  StGlobalDirection BG{BFieldG[0], BFieldG[1], BFieldG[2]};                  PrPP(Make, BG);
  transform_.global_to_local_sector_dir( dirG, segment.dirLS, sector, coorS.row);  PrPP(Make, segment.dirLS);
  transform_.global_to_local_sector_dir(   BG, segment.BLS,   sector, coorS.row);  PrPP(Make, segment.BLS);

  // Distortions
  if (TESTBIT(options_, kDistortion)) {
    float pos[3] = {(float ) coorLT.position.x, (float ) coorLT.position.y, (float ) coorLT.position.z};
    float posMoved[3];
    mag_field_utils_.DoDistortion(pos, posMoved, sector); // input pos[], returns posMoved[]
    coorLT.position = {posMoved[0], posMoved[1], posMoved[2]};       // after distortions
    transform_.local_to_global(coorLT, xyzG);
    PrPP(Make, coorLT);
  }

  transform_.local_to_local_sector(coorLT, segment.coorLS);
  PrPP(Make, segment.coorLS);

  double driftLength = segment.coorLS.position.z + hit.tof * tpcrs::DriftVelocity(sector, cfg_);

  if (driftLength > -1.0 && driftLength <= 0) {
    if ((!IsInner(coorS.row, sector) && driftLength > - cfg_.S<tpcWirePlanes>().outerSectorAnodeWirePadSep) ||
        ( IsInner(coorS.row, sector) && driftLength > - cfg_.S<tpcWirePlanes>().innerSectorAnodeWirePadSep))
      driftLength = std::abs(driftLength);
  }

  segment.coorLS.position.z = driftLength;
  PrPP(Make, segment.coorLS);
  transform_.local_sector_to_hardware(segment.coorLS, segment.Pad, false, false); // don't use T0, don't use Tau
  PrPP(Make, segment.Pad);

  return segment;
}


std::vector<float> Simulator::NumberOfElectronsInCluster(const TF1& heed, float dE, float& dEr)
{
  std::vector<float> rs;

  float dET = dE + dEr;
  dEr = dET;
  float EC;

  while ((EC = const_cast<TF1&>(heed).GetRandom()) < dEr) {
    dEr -= EC;
    rs.push_back(1 - dEr / dET);
  }

  return rs;
}


double Simulator::CalcBaseGain(int sector, int row)
{
  // switch between Inner / Outer Sector paramters
  int iowe = 0;
  if (sector > 12)  iowe += 4;
  if (!IsInner(row, sector)) iowe += 2;

  // Extra correction for simulation with respect to data
  const float* AdditionalMcCorrection = cfg_.S<TpcResponseSimulator>().SecRowCorIW;
  const float* AddSigmaMcCorrection   = cfg_.S<TpcResponseSimulator>().SecRowSigIW;

  double gain = tpcrs::GainCorrection(sector, row, cfg_);
  double gain_x_correctionL = AdditionalMcCorrection[iowe] + row * AdditionalMcCorrection[iowe + 1];
  double gain_x_sigma = AddSigmaMcCorrection[iowe] + row * AddSigmaMcCorrection[iowe + 1];

  gain *= std::exp(-gain_x_correctionL);

  if (gain_x_sigma > 0) gain *= std::exp(gRandom->Gaus(0., gain_x_sigma));

  return gain;
}


double Simulator::CalcLocalGain(int sector, int row, double gain_base, double dedx_corr)
{
  double num_electrons_per_adc = cfg_.S<TpcResponseSimulator>().NoElPerAdc;

  if (num_electrons_per_adc <= 0) {
    if (cfg_.C<St_tpcPadConfigC>().iTPC(sector) && cfg_.C<St_tpcPadConfigC>().IsRowInner(sector, row)) {
      num_electrons_per_adc = cfg_.S<TpcResponseSimulator>().NoElPerAdcX; // iTPC
    }
    else if (cfg_.C<St_tpcPadConfigC>().IsRowInner(sector, row)) {
      num_electrons_per_adc = cfg_.S<TpcResponseSimulator>().NoElPerAdcI; // inner TPX
    }
    else {
      num_electrons_per_adc = cfg_.S<TpcResponseSimulator>().NoElPerAdcO; // outer TPX
    }
  }

  return gain_base / dedx_corr / num_electrons_per_adc;
}


void Simulator::SignalFromSegment(const TrackSegment& segment, TrackHelix track, double gain_local,
  ChargeContainer& binned_charge, int& nP, double& dESum, double& dSSum)
{
  static const double m_e = .51099907e-3;
  static const double eV = 1e-9; // electronvolt in GeV
  static const double cLog10 = std::log(10.);

  double gamma = std::pow(10., segment.tpc_hitC->lgam) + 1;
  double betaGamma = std::sqrt(gamma * gamma - 1.);
  double eKin = -1;
  Coords pxyzG{segment.tpc_hitC->p[0], segment.tpc_hitC->p[1], segment.tpc_hitC->p[2]};
  double bg = segment.mass > 0 ? pxyzG.mag() / segment.mass : 0;

  // special case of stopped electrons
  if (segment.tpc_hitC->particle_id == 3 && segment.tpc_hitC->ds < 0.0050 && segment.tpc_hitC->de < 0) {
    eKin = -segment.tpc_hitC->de;
    gamma = eKin / m_e + 1;
    bg = std::sqrt(gamma * gamma - 1.);
  }

  if (bg > betaGamma) betaGamma = bg;

  gamma = std::sqrt(betaGamma * betaGamma + 1.);
  double Tmax;

  if (segment.mass < 2 * m_e) {
    if (segment.charge > 0) Tmax =       m_e * (gamma - 1);
    else                    Tmax = 0.5 * m_e * (gamma - 1);
  }
  else {
    double r = m_e / segment.mass;
    Tmax = 2 * m_e * betaGamma * betaGamma / (1 + 2 * gamma * r + r * r);
  }

  if (Tmax > cfg_.S<ResponseSimulator>().electron_cutoff_energy)
    Tmax = cfg_.S<ResponseSimulator>().electron_cutoff_energy;

  float dEr = 0;
  double s_low   = -std::abs(segment.tpc_hitC->ds) / 2;
  double s_upper =  std::abs(segment.tpc_hitC->ds) / 2;
  double newPosition = s_low;

  // generate electrons: No. of primary clusters per cm
  double NP = GetNoPrimaryClusters(betaGamma, segment.charge); // per cm

  do {// Clusters
    float dS = 0;
    float dE = 0;

    if (eKin >= 0.0) {
      if (eKin == 0.0) break;

      double gamma = eKin / m_e + 1;
      double bg = std::sqrt(gamma * gamma - 1.);
      Tmax = 0.5 * m_e * (gamma - 1);

      if (Tmax <= cfg_.S<TpcResponseSimulator>().W / 2 * eV) break;

      NP = GetNoPrimaryClusters(betaGamma, segment.charge);
      dE = std::exp(cLog10 * mdNdEL10->GetRandom());
    }
    else {
      if (segment.charge) {
        dS = - std::log(gRandom->Rndm()) / NP;

        if (mdNdEL10) dE = std::exp(cLog10 * mdNdEL10->GetRandom());
        else          dE = cfg_.S<TpcResponseSimulator>().W *
                           gRandom->Poisson(cfg_.S<TpcResponseSimulator>().Cluster);
      }
      else { // charge == 0 geantino
        // for LASERINO assume dE/dx = 25 keV/cm;
        dE = 10; // eV
        dS = dE * eV / (std::abs(segment.tpc_hitC->de / segment.tpc_hitC->ds));
      }
    }

    double E = dE * eV;
    newPosition += dS;

    if (newPosition > s_upper) break;

    if (dE < cfg_.S<TpcResponseSimulator>().W / 2 || E > Tmax) continue;

    if (eKin > 0) {
      if (eKin >= E) {eKin -= E;}
      else {E = eKin; eKin = 0; dE = E / eV;}
    }

    dESum += dE;
    dSSum += dS;
    nP++;
    double xRange = ElectronRange(dE, dEr);

    std::vector<float> rs = NumberOfElectronsInCluster(mHeed, dE, dEr);

    if (!rs.size()) continue;

    Coords xyzC = track.at(newPosition);

    LoopOverElectronsInCluster(rs, segment, binned_charge, xRange, xyzC, gain_local);
  }
  while (true);   // Clusters
}


void Simulator::LoopOverElectronsInCluster(
  std::vector<float> rs, const TrackSegment &segment, ChargeContainer& binned_charge,
  double xRange, Coords xyzC, double gain_local)
{
  int sector = segment.Pad.sector;
  int row    = segment.Pad.row;
  double OmegaTau = cfg_.S<TpcResponseSimulator>().OmegaTau *
                      segment.BLS.position.z / 5.0; // from diffusion 586 um / 106 um at B = 0/ 5kG
  double driftLength = std::abs(segment.coorLS.position.z);
  double D = 1. + OmegaTau * OmegaTau;
  double SigmaL = cfg_.S<TpcResponseSimulator>().longitudinalDiffusion * std::sqrt(driftLength);
  double SigmaT = cfg_.S<TpcResponseSimulator>().transverseDiffusion * std::sqrt(driftLength / D);
  double sigmaJitterX = IsInner(row, sector) ?  cfg_.S<TpcResponseSimulator>().SigmaJitterXI :
                                                cfg_.S<TpcResponseSimulator>().SigmaJitterXO;
  if (sigmaJitterX > 0) {
    SigmaT = std::sqrt(SigmaT * SigmaT + sigmaJitterX * sigmaJitterX);
  }

  // Dummy call to keep the same random number sequence
  gRandom->Rndm();
  double rX, rY;
  int WireIndex = 0;

  InOut io = IsInner(row, sector) ? kInner : kOuter;

  Coords unit = segment.dirLS.position.unit();
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
      for (int i=0; i<3; i++) xyzE.position[i] += xyzR[i];
    }

    double y = xyzE.position.y;
    double alphaVariation = InnerAlphaVariation[sector - 1];

    // Transport to wire
    if (y <= lastInnerSectorAnodeWire) {
      WireIndex = tpcrs::irint((y - firstInnerSectorAnodeWire) / anodeWirePitch) + 1;
      if (cfg_.C<St_tpcPadConfigC>().iTPC(sector)) {// two first and two last wires are removed, and 3rd wire is fat wiere
        if (WireIndex <= 3 || WireIndex >= numberOfInnerSectorAnodeWires - 3) continue;
      }
      else {   // old TPC the first and last wires are fat ones
        if (WireIndex <= 1 || WireIndex >= numberOfInnerSectorAnodeWires) continue;
      }
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
    double tanLorentz = OmegaTau / cfg_.S<TpcResponseSimulator>().OmegaTauScaleO;

    if (y < firstOuterSectorAnodeWire) tanLorentz = OmegaTau / cfg_.S<TpcResponseSimulator>().OmegaTauScaleI;

    xOnWire += distFocused * tanLorentz; // tanLorentz near wires taken from comparison with data
    zOnWire += std::abs(distFocused);

    if (! iGroundWire ) gain_gas *= std::exp( alphaVariation);
    else                gain_gas *= std::exp(-alphaVariation);

    double dY     = mChargeFraction[io][sector - 1].GetXmax();
    double yLmin  = yOnWire - dY;
    double yLmax  = yOnWire + dY;

    int    rowMin = transform_.rowFromLocalY(yLmin, sector);
    int    rowMax = transform_.rowFromLocalY(yLmax, sector);
    double yRmin  = transform_.yFromRow(sector, rowMin) - cfg_.C<St_tpcPadConfigC>().PadLengthAtRow(sector, rowMin) / 2;
    double yRmax  = transform_.yFromRow(sector, rowMax) + cfg_.C<St_tpcPadConfigC>().PadLengthAtRow(sector, rowMax) / 2;

    if (yRmin > yLmax || yRmax < yLmin) {
      continue;
    }

    GenerateSignal(segment, rowMin, rowMax,
                   &mShaperResponses[io][sector - 1], binned_charge, gain_local * gain_gas);
  }  // electrons in Cluster
}


void Simulator::GenerateSignal(const TrackSegment &segment, int rowMin, int rowMax,
  TF1F* shaper, ChargeContainer& binned_charge, double gain_local_gas)
{
  int sector = segment.Pad.sector;
  int row    = segment.Pad.row;
  double sigmaJitterT = (IsInner(row, sector) ? cfg_.S<TpcResponseSimulator>().SigmaJitterTI :
                                                cfg_.S<TpcResponseSimulator>().SigmaJitterTO);

  for (int row = rowMin; row <= rowMax; row++) {
    if (cfg_.C<St_tpcPadConfigC>().numberOfRows(sector) == 45) { // ! iTpx
      if ( !cfg_.C<St_tpcRDOMasksC>().isRowOn(sector, row) ) continue;
      if ( !cfg_.C<St_tpcAnodeHVavgC>().livePadrow(sector, row) )  continue;
    }

    StTpcLocalSectorCoordinate xyzW{xOnWire, yOnWire, zOnWire, sector, row};
    StTpcPadCoordinate Pad;
    transform_.local_sector_to_hardware(xyzW, Pad, false, false); // don't use T0, don't use Tau
    float bin = Pad.timeBucket;//L  - 1; // K
    int binT = tpcrs::irint(bin); //L bin;//K tpcrs::irint(bin);// J bin; // I tpcrs::irint(bin);

    if (binT < 0 || binT >= max_timebins_) continue;

    double dT = bin - binT + cfg_.S<TpcResponseSimulator>().T0offset;
    dT += IsInner(row, sector) ? cfg_.S<TpcResponseSimulator>().T0offsetI :
                                 cfg_.S<TpcResponseSimulator>().T0offsetO;

    if (sigmaJitterT) dT += gRandom->Gaus(0, sigmaJitterT);

    InOut io = IsInner(row, sector) ? kInner : kOuter;

    double delta_y = transform_.yFromRow(sector, row) - yOnWire;
    double YDirectionCoupling = mChargeFraction[io][sector - 1].GetSaveL(delta_y);

    if (YDirectionCoupling < cfg_.S<ResponseSimulator>().min_signal) continue;

    float padX = Pad.pad;
    int CentralPad = tpcrs::irint(padX);
    int PadsAtRow = cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(sector, row);

    if (CentralPad < 1 || CentralPad > PadsAtRow) continue;

    int DeltaPad = tpcrs::irint(mPadResponseFunction[io][sector - 1].GetXmax()) + 1;
    int padMin   = std::max(CentralPad - DeltaPad, 1);
    int padMax   = std::min(CentralPad + DeltaPad, PadsAtRow);
    int Npads    = std::min(padMax - padMin + 1, static_cast<int>(kPadMax));
    double xPadMin = padMin - padX;

    static double XDirectionCouplings[kPadMax];
    mPadResponseFunction[io][sector - 1].GetSaveL(Npads, xPadMin, XDirectionCouplings);

    for (int pad = padMin; pad <= padMax; pad++) {
      double gain = gain_local_gas;
      double dt = dT;

      if ( !TESTBIT(options_, kGAINOAtALL) ) {
        gain *= cfg_.C<St_tpcPadGainT0BC>().Gain(sector, row, pad);

        if (gain <= 0.0) continue;

        dt -= cfg_.C<St_tpcPadGainT0BC>().T0(sector, row, pad);
      }

      double XYcoupling = gain * XDirectionCouplings[pad - padMin] * YDirectionCoupling;

      if (XYcoupling < cfg_.S<ResponseSimulator>().min_signal) continue;

      int tbin_first = std::max(0, binT + tpcrs::irint(dt + shaper->GetXmin() - 0.5));
      int tbin_last  = std::min(max_timebins_ - 1, binT + tpcrs::irint(dt + shaper->GetXmax() + 0.5));
      int num_tbins  = std::min(tbin_last - tbin_first + 1, static_cast<int>(kTimeBacketMax));

      static double TimeCouplings[kTimeBacketMax];
      shaper->GetSaveL(num_tbins, tbin_first - binT - dt, TimeCouplings);

      int index = max_timebins_ * ((row - 1) * max_pads_ + pad - 1) + tbin_first;

      for (int itbin = tbin_first; itbin <= tbin_last; itbin++, index++) {
        double signal = XYcoupling * TimeCouplings[itbin - tbin_first];

        if (signal < cfg_.S<ResponseSimulator>().min_signal)  continue;

        binned_charge[index].Sum += signal;

        // Record truth ID of the MC particle produced this signal
        if ( segment.tpc_hitC->track_id ) {
          if ( !binned_charge[index].TrackId )
            binned_charge[index].TrackId = segment.tpc_hitC->track_id;
          else { // switch TrackId, works only for 2 tracks, more tracks ?
            if ( binned_charge[index].TrackId != segment.tpc_hitC->track_id && binned_charge[index].Sum < 2 * signal)
              binned_charge[index].TrackId = segment.tpc_hitC->track_id;
          }
        }
      } // time
    } // pad limits
  } // row limits
}


double Simulator::dEdxCorrection(const TrackSegment &segment)
{
  dEdxY2_t CdEdx;
  memset(&CdEdx, 0, sizeof(dEdxY2_t));
  CdEdx.DeltaZ  = 5.2;
  CdEdx.QRatio  = -2;
  CdEdx.QRatioA = -2.;
  CdEdx.QSumA   = 0;
  CdEdx.sector  = segment.Pad.sector;
  CdEdx.row     = segment.Pad.row;
  CdEdx.pad     = tpcrs::irint(segment.Pad.pad);
  CdEdx.edge    = CdEdx.pad;

  if (CdEdx.edge > 0.5 * cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(CdEdx.sector, CdEdx.row))
    CdEdx.edge += 1 - cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(CdEdx.sector, CdEdx.row);

  CdEdx.F.dE   = 1;
  CdEdx.F.dx   = std::abs(segment.tpc_hitC->ds);
  CdEdx.xyz[0] = segment.coorLS.position.x;
  CdEdx.xyz[1] = segment.coorLS.position.y;
  CdEdx.xyz[2] = segment.coorLS.position.z;
  double probablePad = cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(CdEdx.sector, CdEdx.row) / 2;
  double pitch = IsInner(CdEdx.row, CdEdx.sector) ? cfg_.C<St_tpcPadConfigC>().innerSectorPadPitch(CdEdx.sector) :
                                                    cfg_.C<St_tpcPadConfigC>().outerSectorPadPitch(CdEdx.sector);
  double PhiMax = std::atan2(probablePad * pitch, cfg_.C<St_tpcPadConfigC>().radialDistanceAtRow(CdEdx.sector, CdEdx.row));
  CdEdx.PhiR    = std::atan2(CdEdx.xyz[0], CdEdx.xyz[1]) / PhiMax;
  CdEdx.xyzD[0] = segment.dirLS.position.x;
  CdEdx.xyzD[1] = segment.dirLS.position.y;
  CdEdx.xyzD[2] = segment.dirLS.position.z;
  CdEdx.zG      = CdEdx.xyz[2];
  CdEdx.Zdc     = cfg_.S<trigDetSums>().zdcX;
  CdEdx.ZdriftDistance = segment.coorLS.position.z; // drift length
  CdEdx.ZdriftDistanceO2 = CdEdx.ZdriftDistance * cfg_.S<tpcGas>().ppmOxygenIn;

  return dEdx_correction_.dEdxCorrection(CdEdx) ? 1 : CdEdx.F.dE;
}

#undef PrPP
