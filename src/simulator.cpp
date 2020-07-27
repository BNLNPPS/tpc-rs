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
  dEdx_correction_(cfg, dEdxCorr::kAll & ~dEdxCorr::kAdcCorrection & ~dEdxCorr::kAdcCorrectionMDF & ~dEdxCorr::kdXCorrection, 0),
  digi_(cfg_),
  alpha_gain_variations_()
{
  //SETBIT(options_, kHEED);
  SETBIT(options_, kBICHSEL); // Default is Bichsel
  SETBIT(options_, kdEdxCorr);
  //SETBIT(options_, kDistortion);

  if (TESTBIT(options_, kBICHSEL)) {
    LOG_INFO << "Model for dE/dx simulation: Bichsel\n";

    TFile inner(cfg_.Locate("dNdE_Bichsel.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(cfg_.Locate("dNdx_Bichsel.root").c_str());
    mdNdx = (TH1D*) outer.Get("dNdx"); assert(mdNdx);
    mdNdx->SetDirectory(0);
  }
  else if (TESTBIT(options_, kHEED)) {
    LOG_INFO << "Model for dE/dx simulation: Heed\n";

    TFile inner(cfg_.Locate("dNdx_Heed.root").c_str());
    mdNdEL10 = (TH1D*) inner.Get("dNdEL10"); assert(mdNdEL10);
    mdNdEL10->SetDirectory(0);

    TFile outer(cfg_.Locate("dNdx_Heed.root").c_str());
    mdNdxL10 = (TH1D*) outer.Get("dNdxL10"); assert(mdNdxL10);
    mdNdxL10->SetDirectory(0);
  }
  else {
    LOG_INFO << "Error: Model for dE/dx simulation not set\n";
  }

  double t0IO[2];
  InitAlphaGainVariations(t0IO);

  double TimeBinWidth = 1. / cfg_.S<starClockOnl>().frequency;

  // Shapers
  double timeBinMin = -0.5;
  double timeBinMax = 44.5;

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
        io == kInner ? cfg_.S<tpcPadPlanes>().innerSectorPadPitch :
                       cfg_.S<tpcPadPlanes>().outerSectorPadPitch
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


void Simulator::InitAlphaGainVariations(double t0IO[2])
{
  alpha_gain_variations_.resize(num_sectors_*2, 0);

  double anode_wire_pitch  = cfg_.S<tpcWirePlanes>().anodeWirePitch;
  double anode_wire_radius = cfg_.S<tpcWirePlanes>().anodeWireRadius;
  double cathod_anode_gap[2] = {0.2, 0.4};
  std::vector<double> avg_anode_voltage(num_sectors_*2, 0);

  for (int sector = 1; sector <= num_sectors_; sector++)
  {
    int n_inner = 0, n_outer = 0;
    for (int row = 1; row <= cfg_.S<tpcPadPlanes>().padRows; row++) {
      if (tpcrs::IsInner(row, cfg_)) {
        n_inner++;
        avg_anode_voltage[num_sectors_*0 + sector - 1] += cfg_.C<St_tpcAnodeHVavgC>().voltagePadrow(sector, row);
      }
      else {
        n_outer++;
        avg_anode_voltage[num_sectors_*1 + sector - 1] += cfg_.C<St_tpcAnodeHVavgC>().voltagePadrow(sector, row);
      }
    }

    avg_anode_voltage[num_sectors_*0 + sector - 1] /= n_inner;
    avg_anode_voltage[num_sectors_*1 + sector - 1] /= n_outer;

    for (int io = kInner; io <= kOuter; io++)
    {
      double volt   = avg_anode_voltage[num_sectors_*io + sector - 1];
      double charge = InducedCharge(anode_wire_pitch, cathod_anode_gap[io],
                                    anode_wire_radius, volt, t0IO[io]);

      alpha_gain_variations_[num_sectors_*io + sector - 1] = charge;
    }
  }
}


template<>
void Simulator::Simulate(SimuHitIt first_hit, SimuHitIt last_hit, DigiInserter digi_data, DistInserter distorted)
{
  static int nCalls = 0;
  gRandom->SetSeed(2345 + nCalls++);

  double vminI = cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min;
  double vminO = cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min;

  cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min = -500;
  cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min = -500;

  LOG_INFO << "Reset min for gain Correction to I/O\t"
           << cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min
           << "\t"
           << cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min
           << " (V)\n";

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

  unsigned int sector = 1;
  for (auto segments_in_sector : segments_by_sector) {
    int nHitsInTheSector = 0;
    ChargeContainer binned_charge(digi_.total_timebins(), {0, 0});

    for (TrackSegment& segment : segments_in_sector) {
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
        segment.coorLS.position = {track.at(sR).x, track.at(sR).y, track.at(sR).z};
        transform_.local_sector_to_hardware(segment.coorLS, segment.Pad, false, false); // don't use T0, don't use Tau
      }

      int nP = 0;
      double dESum = 0;
      double dSSum = 0;

      SignalFromSegment(segment, track, gain_local, binned_charge, nP, dESum, dSSum);

      *distorted = tpcrs::DistortedHit{
        {segment.coorLS.position.x, segment.coorLS.position.y, segment.coorLS.position.z},
        {segment.dirLS.position.x,  segment.dirLS.position.y,  segment.dirLS.position.z},
        dESum * 1e-9, // electronvolt in GeV
        dSSum,
        nP
      };

      nHitsInTheSector++;
    } // end do loop over segments for a given particle

    if (nHitsInTheSector) {
      DigitizeSector(sector, binned_charge, digi_data);
    }

    sector++;
  }

  cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min = vminI;
  cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min = vminO;

  LOG_INFO << "Reset min for gain Correction to I/O\t"
           << cfg_.C<St_tpcGainCorrectionC>().Struct(1)->min
           << "\t"
           << cfg_.C<St_tpcGainCorrectionC>().Struct(0)->min
           << " (V)\n";
}


void Simulator::CreateTrackSegments(SimuHitIt first_hit, SimuHitIt last_hit, std::vector<TrackSegment>& segments)
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


void Simulator::DigitizeSector(unsigned int sector, const ChargeContainer& binned_charge, DigiInserter digi_data)
{
  double pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX;
  double ped = cfg_.S<TpcResponseSimulator>().AveragePedestal;

  std::vector<short> ADCs_(binned_charge.size(), 0);

  auto bc = binned_charge.begin();
  auto adcs_iter = ADCs_.begin();

  for (auto ch = digi_.channels.begin(); ch != digi_.channels.end(); ch += max_timebins_)
  {
    double gain = cfg_.S<tpcPadGainT0>().Gain[sector-1][ch->row-1][ch->pad-1];

    if (gain <= 0) {
      bc        += max_timebins_;
      adcs_iter += max_timebins_;
      continue;
    }

    for (int i=0; i != max_timebins_; ++i, ++bc, ++adcs_iter)
    {
      int adc = int(bc->Sum / gain + gRandom->Gaus(ped, pedRMS) - ped);
      // Zero negative values
      adc = adc & ~(adc >> 31);
      // Select minimum between adc and 1023, i.e. overflow at 1023
      adc = adc - !(((adc - 1023) >> 31) & 0x1) * (adc - 1023);

      *adcs_iter = adc;
    }
  }

  for (auto adcs_iter = ADCs_.begin(); adcs_iter != ADCs_.end(); adcs_iter += max_timebins_)
  {
    SimulateAltro(adcs_iter, adcs_iter + max_timebins_, true);
  }

  auto ch = digi_.channels.begin();
  adcs_iter = ADCs_.begin();

  for (auto bc = binned_charge.begin(); bc != binned_charge.end(); ++bc, ++ch, ++adcs_iter)
  {
    if (*adcs_iter == 0) continue;
    *digi_data = tpcrs::DigiHit{sector, ch->row, ch->pad, ch->timebin, *adcs_iter, bc->TrackId};
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
  const double B  = 30e-3; // 1/V
  const double E0 = 20e3; // V/cm
  const double mu = 2.26; // cm**2/V/sec CH4+ mobility
  // const double mu = 1.87; // cm**2/V/sec Ar+ mobility
  double alpha[2] = {-26., -70.};
  // E.Mathieson (3.2b), V.Chernyatin said that it should be used this (Weber ) approximation 07/09/08
  double rc = s / (2 * M_PI) * std::exp(M_PI * h / s);
  double C  = 1. / (2 * std::log(rc / ra));
  double E  = 2 * M_PI * C * Va / s;
  // Gain variation: M = M0*(1 - k*cos(2*alpha))
  double k = 2 * B / 3.*std::pow((M_PI / E0 / s), 2) * std::pow(C * Va, 3);
  // Induced charge variation
  t0 = ra * ra / (4 * mu * C * Va);
  double Tav = t0 * h / s / (2 * M_PI * C);
  double t = 180e-9;
  double rp = std::sqrt(1. + t / t0);
  double Aconstant = rp * ra / (2 * h);
  double Bconstant = C / 2 * std::log(1 + t / t0);
  double Gains[2];

  for (int i = 0; i < 2; i++) {
    Gains[i] = Aconstant * std::sin(M_PI / 180 * alpha[i]) + Bconstant;
  }

  double GainsAv = std::sqrt(Gains[0] * Gains[1]);
  double r = 0;

  for (int i = 0; i < 2; i++) {
    r = std::log(Gains[i] / GainsAv);
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
  StGlobalCoordinate xyzG{hit.x[0], hit.x[1], hit.x[2]};
  segment.simu_hit = &hit;
  ParticleProperties(hit.particle_id, segment.charge, segment.mass);

  StTpcLocalSectorCoordinate coorS;
  // GlobalCoord -> LocalSectorCoord. This transformation can result in a row
  // that is not the same as (volId % 100)
  transform_.global_to_local_sector(xyzG, coorS, sector, 0);
  StTpcLocalCoordinate coorLT;  // before distortions
  transform_.global_to_local(xyzG, coorLT, sector, coorS.row);

  // move up, calculate field at center of TPC
  static float BFieldG[3];
  mag_field_utils_.BFieldTpc(hit.x, BFieldG);
  // distortion and misalignment
  // replace pxy => direction and try linear extrapolation
  Coords pxyzG{hit.p[0], hit.p[1], hit.p[2]};
  StGlobalDirection dirG{pxyzG.unit()};
  StGlobalDirection BG{BFieldG[0], BFieldG[1], BFieldG[2]};
  transform_.global_to_local_sector_dir( dirG, segment.dirLS, sector, coorS.row);
  transform_.global_to_local_sector_dir(   BG, segment.BLS,   sector, coorS.row);

  // Distortions
  if (TESTBIT(options_, kDistortion)) {
    float pos[3] = {(float ) coorLT.position.x, (float ) coorLT.position.y, (float ) coorLT.position.z};
    float posMoved[3];
    mag_field_utils_.DoDistortion(pos, posMoved, sector); // input pos[], returns posMoved[]
    coorLT.position = {posMoved[0], posMoved[1], posMoved[2]};       // after distortions
    transform_.local_to_global(coorLT, xyzG);
  }

  transform_.local_to_local_sector(coorLT, segment.coorLS);

  double driftLength = segment.coorLS.position.z + hit.tof * tpcrs::DriftVelocity(sector, cfg_);

  if (driftLength > -1.0 && driftLength <= 0) {
    if ((!tpcrs::IsInner(coorS.row, cfg_) && driftLength > - cfg_.S<tpcWirePlanes>().outerSectorAnodeWirePadSep) ||
        ( tpcrs::IsInner(coorS.row, cfg_) && driftLength > - cfg_.S<tpcWirePlanes>().innerSectorAnodeWirePadSep))
      driftLength = std::abs(driftLength);
  }

  segment.coorLS.position.z = driftLength;
  transform_.local_sector_to_hardware(segment.coorLS, segment.Pad, false, false); // don't use T0, don't use Tau

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
  if (!tpcrs::IsInner(row, cfg_)) iowe += 2;

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
  return gain_base / dedx_corr / cfg_.S<TpcResponseSimulator>().NoElPerAdc;
}


void Simulator::SignalFromSegment(const TrackSegment& segment, TrackHelix track, double gain_local,
  ChargeContainer& binned_charge, int& nP, double& dESum, double& dSSum)
{
  static const double m_e = .51099907e-3;
  static const double eV = 1e-9; // electronvolt in GeV
  static const double cLog10 = std::log(10.);

  double gamma = std::pow(10., segment.simu_hit->lgam) + 1;
  double betaGamma = std::sqrt(gamma * gamma - 1.);
  double eKin = -1;
  Coords pxyzG{segment.simu_hit->p[0], segment.simu_hit->p[1], segment.simu_hit->p[2]};
  double bg = segment.mass > 0 ? pxyzG.mag() / segment.mass : 0;

  // special case of stopped electrons
  if (segment.simu_hit->particle_id == 3 && segment.simu_hit->ds < 0.0050 && segment.simu_hit->de < 0) {
    eKin = -segment.simu_hit->de;
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
  double s_low   = -std::abs(segment.simu_hit->ds) / 2;
  double s_upper =  std::abs(segment.simu_hit->ds) / 2;
  double newPosition = s_low;

  // generate electrons: No. of primary clusters per cm
  double NP = GetNoPrimaryClusters(betaGamma, segment.charge); // per cm

  do {// Clusters
    float dS = 0;
    float dE = 0;

    if (eKin >= 0.0) {
      if (eKin == 0.0) break;

      double gamma = eKin / m_e + 1;
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
        dS = dE * eV / (std::abs(segment.simu_hit->de / segment.simu_hit->ds));
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
  double sigmaJitterX = tpcrs::IsInner(row, cfg_) ?  cfg_.S<TpcResponseSimulator>().SigmaJitterXI :
                                                cfg_.S<TpcResponseSimulator>().SigmaJitterXO;
  if (sigmaJitterX > 0) {
    SigmaT = std::sqrt(SigmaT * SigmaT + sigmaJitterX * sigmaJitterX);
  }

  // Dummy call to keep the same random number sequence
  gRandom->Rndm();
  double rX, rY;

  InOut io = tpcrs::IsInner(row, cfg_) ? kInner : kOuter;

  Coords unit = segment.dirLS.position.unit();
  double L2L[9] = {unit.z,                  - unit.x*unit.z, unit.x,
                   unit.x,                  - unit.y*unit.z, unit.y,
                   0.0,       unit.x*unit.x + unit.y*unit.y, unit.z};

  for (int ie = 0; ie < rs.size(); ie++)
  {
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

    bool missed_readout = false;
    bool is_ground_wire = false;

    Coords at_readout = TransportToReadout(xyzE.position, OmegaTau, missed_readout, is_ground_wire);

    if (missed_readout) continue;

    double alphaVariation = (xyzE.position.y <= cfg_.S<tpcWirePlanes>().lastInnerSectorAnodeWire) ?
                             alpha_gain_variations_[num_sectors_*0 + sector - 1] :
                             alpha_gain_variations_[num_sectors_*1 + sector - 1];

    if (!is_ground_wire) gain_gas *= std::exp( alphaVariation);
    else                 gain_gas *= std::exp(-alphaVariation);

    double dY     = mChargeFraction[io][sector - 1].GetXmax();
    double yLmin  = at_readout.y - dY;
    double yLmax  = at_readout.y + dY;

    int    rowMin = transform_.rowFromLocalY(yLmin, sector);
    int    rowMax = transform_.rowFromLocalY(yLmax, sector);
    double yRmin  = transform_.yFromRow(sector, rowMin) - cfg_.C<St_tpcPadConfigC>().PadLengthAtRow(sector, rowMin) / 2;
    double yRmax  = transform_.yFromRow(sector, rowMax) + cfg_.C<St_tpcPadConfigC>().PadLengthAtRow(sector, rowMax) / 2;

    if (yRmin > yLmax || yRmax < yLmin) {
      continue;
    }

    GenerateSignal(segment, at_readout, rowMin, rowMax,
                   &mShaperResponses[io][sector - 1], binned_charge, gain_local * gain_gas);
  }  // electrons in Cluster
}


void Simulator::GenerateSignal(const TrackSegment &segment, Coords at_readout, int rowMin, int rowMax,
  TF1F* shaper, ChargeContainer& binned_charge, double gain_local_gas)
{
  int sector = segment.Pad.sector;
  int row    = segment.Pad.row;
  double sigmaJitterT = (tpcrs::IsInner(row, cfg_) ? cfg_.S<TpcResponseSimulator>().SigmaJitterTI :
                                                cfg_.S<TpcResponseSimulator>().SigmaJitterTO);

  for (unsigned row = rowMin; row <= rowMax; row++) {
    if ( !cfg_.C<St_tpcAnodeHVavgC>().livePadrow(sector, row) )  continue;

    StTpcLocalSectorCoordinate xyzW{at_readout.x, at_readout.y, at_readout.z, sector, row};
    StTpcPadCoordinate Pad;
    transform_.local_sector_to_hardware(xyzW, Pad, false, false); // don't use T0, don't use Tau
    float bin = Pad.timeBucket;//L  - 1; // K
    int binT = tpcrs::irint(bin); //L bin;//K tpcrs::irint(bin);// J bin; // I tpcrs::irint(bin);

    if (binT < 0 || binT >= max_timebins_) continue;

    double dT = bin - binT + cfg_.S<TpcResponseSimulator>().T0offset;
    dT += tpcrs::IsInner(row, cfg_) ? cfg_.S<TpcResponseSimulator>().T0offsetI :
                                 cfg_.S<TpcResponseSimulator>().T0offsetO;

    if (sigmaJitterT) dT += gRandom->Gaus(0, sigmaJitterT);

    InOut io = tpcrs::IsInner(row, cfg_) ? kInner : kOuter;

    double delta_y = transform_.yFromRow(sector, row) - at_readout.y;
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

    for (unsigned pad = padMin; pad <= padMax; pad++) {
      double gain = gain_local_gas;
      double dt = dT;

      if ( !TESTBIT(options_, kGAINOAtALL) ) {
        gain *= cfg_.S<tpcPadGainT0>().Gain[sector-1][row-1][pad-1];

        if (gain <= 0.0) continue;

        dt -= cfg_.S<tpcPadGainT0>().T0[sector-1][row-1][pad-1];
      }

      double XYcoupling = gain * XDirectionCouplings[pad - padMin] * YDirectionCoupling;

      if (XYcoupling < cfg_.S<ResponseSimulator>().min_signal) continue;

      int tbin_first = std::max(0, binT + tpcrs::irint(dt + shaper->GetXmin() - 0.5));
      int tbin_last  = std::min(max_timebins_ - 1, binT + tpcrs::irint(dt + shaper->GetXmax() + 0.5));
      int num_tbins  = std::min(tbin_last - tbin_first + 1, static_cast<int>(kTimeBacketMax));

      static double TimeCouplings[kTimeBacketMax];
      shaper->GetSaveL(num_tbins, tbin_first - binT - dt, TimeCouplings);

      int index = max_timebins_ * (digi_.total_pads(row) + pad - 1) + tbin_first;

      for (unsigned itbin = tbin_first; itbin <= tbin_last; itbin++, index++) {
        double signal = XYcoupling * TimeCouplings[itbin - tbin_first];

        if (signal < cfg_.S<ResponseSimulator>().min_signal) continue;

        binned_charge[index] += {signal, segment.simu_hit->track_id};
      } // time
    } // pad limits
  } // row limits
}


Coords Simulator::TransportToReadout(const Coords c, double OmegaTau, bool& missed_readout, bool& is_ground_wire)
{
  Coords readout;
  missed_readout = false;

  // Transport to wire
  double firstOuterSectorAnodeWire = cfg_.S<tpcWirePlanes>().firstOuterSectorAnodeWire;
  double anodeWirePitch            = cfg_.S<tpcWirePlanes>().anodeWirePitch;
  int wire_index = 0;

  if (c.y <= cfg_.S<tpcWirePlanes>().lastInnerSectorAnodeWire) {
    double firstInnerSectorAnodeWire = cfg_.S<tpcWirePlanes>().firstInnerSectorAnodeWire;
    wire_index = tpcrs::irint((c.y - firstInnerSectorAnodeWire) / anodeWirePitch) + 1;
    // In TPC the first and last wires are fat ones
    if (wire_index <= 1 || wire_index >= cfg_.S<tpcWirePlanes>().numInnerSectorAnodeWires)
      missed_readout = true;
    readout.y = firstInnerSectorAnodeWire + (wire_index - 1) * anodeWirePitch;
  }
  else {
    wire_index = tpcrs::irint((c.y - firstOuterSectorAnodeWire) / anodeWirePitch) + 1;
    // In TPC the first and last wires are fat ones
    if (wire_index <= 1 || wire_index >= cfg_.S<tpcWirePlanes>().numOuterSectorAnodeWires)
      missed_readout = true;
    readout.y = firstOuterSectorAnodeWire + (wire_index - 1) * anodeWirePitch;
  }

  double distance_to_wire = c.y - readout.y; // Calculated effective distance to wire affected by Lorentz shift
  // Grid plane (1 mm spacing) focusing effect + Lorentz angle in drift volume
  int iGridWire = int(std::abs(10.*distance_to_wire));
  double dist2Grid = std::copysign(0.05 + 0.1 * iGridWire, distance_to_wire); // [cm]
  // Ground plane (1 mm spacing) focusing effect
  int iGroundWire = int(std::abs(10.*dist2Grid));
  double distFocused = std::copysign(0.05 + 0.1 * iGroundWire, dist2Grid);

  // OmegaTau near wires taken from comparison with data
  double tanLorentz = (c.y < firstOuterSectorAnodeWire) ? OmegaTau / cfg_.S<TpcResponseSimulator>().OmegaTauScaleI :
                                                          OmegaTau / cfg_.S<TpcResponseSimulator>().OmegaTauScaleO;

  readout.x = c.x + distFocused * tanLorentz; // tanLorentz near wires taken from comparison with data
  readout.z = c.z + std::abs(distFocused);

  is_ground_wire = iGroundWire != 0;

  return readout;
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
  CdEdx.F.dx   = std::abs(segment.simu_hit->ds);
  CdEdx.xyz[0] = segment.coorLS.position.x;
  CdEdx.xyz[1] = segment.coorLS.position.y;
  CdEdx.xyz[2] = segment.coorLS.position.z;
  double probablePad = cfg_.C<St_tpcPadConfigC>().numberOfPadsAtRow(CdEdx.sector, CdEdx.row) / 2;
  double pitch = tpcrs::IsInner(CdEdx.row, cfg_) ? cfg_.S<tpcPadPlanes>().innerSectorPadPitch :
                                                   cfg_.S<tpcPadPlanes>().outerSectorPadPitch;
  double PhiMax = std::atan2(probablePad * pitch, tpcrs::RadialDistanceAtRow(CdEdx.row, cfg_));
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
