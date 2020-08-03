#include "bichsel.h"
#include "dedx_correction.h"
#include "logger.h"
#include "struct_containers.h"


StTpcdEdxCorrection::StTpcdEdxCorrection(const tpcrs::Configurator& cfg, int options, int debug) :
  cfg_(cfg),
  options_(options),
  corrections_{
    {"UnCorrected",            ""									, 0, 0, 0},
    {"TpcAdcCorrectionB",      "ADC/Clustering nonlinearity correction"					, &cfg_.C<St_TpcAdcCorrectionBC>(), 0, 0},
    {"TpcEdge",                "Gain on distance from Chamber edge"					, &cfg_.C<St_TpcEdgeC>(), 0, 0},
    {"TpcAdcCorrectionMDF",    "ADC/Clustering nonlinearity correction MDF"				, &cfg_.C<St_TpcAdcCorrectionMDF>(), 0, 0},
    {"TpcdCharge",             "ADC/Clustering undershoot correction"					, &cfg_.C<St_TpcdChargeC>(), 0, 0},
    {"TpcrCharge",             "ADC/Clustering rounding correction"					, &cfg_.C<St_TpcrChargeC>(), 0, 0},
    {"TpcCurrentCorrection",   "Correction due to sagg of Voltage due to anode current"			, &cfg_.C<St_TpcCurrentCorrectionC>(), 0, 0},
    {"TpcSecRowB",             "Gas gain correction for sector/row"					, &cfg_.C<St_TpcSecRowBC>(), 0, 0},
    {"TpcSecRowC",             "Additional Gas gain correction for sector/row"				, &cfg_.C<St_TpcSecRowCC>(), 0, 0},
    {"TpcRowQ",                "Gas gain correction for row versus accumulated charge,"			, &cfg_.C<St_TpcRowQC>(), 0, 0},
    {"tpcPressureB",           "Gain on Gas Density due to Pressure"					, &cfg_.C<St_tpcPressureBC>(), 0, 0},
    {"tpcTime"       ,         "Unregognized time dependce"						, &cfg_.C<St_tpcTimeDependenceC>(), 0, 0},
    {"TpcDriftDistOxygen",     "Correction for Electron Attachment due to O2"				, &cfg_.C<St_TpcDriftDistOxygenC>(), 0, 0},
    {"TpcMultiplicity",        "Global track multiplicity dependence"					, &cfg_.C<St_TpcMultiplicityC>(), 0, 0},
    {"TpcZCorrectionB",        "Variation on drift distance"						, &cfg_.C<St_TpcZCorrectionBC>(), 0, 0},
    {"tpcMethaneIn",           "Gain on Methane content"						, &cfg_.C<St_tpcMethaneInC>(), 0, 0},
    {"tpcGasTemperature",      "Gain on Gas Dens. due to Temperature"					, &cfg_.C<St_tpcGasTemperatureC>(), 0, 0},
    {"tpcWaterOut",            "Gain on Water content"							, &cfg_.C<St_tpcWaterOutC>(), 0, 0},
    {"TpcSpaceCharge",         "Gain on space charge near the wire"					, &cfg_.C<St_TpcSpaceChargeC>(), 0, 0},
    {"TpcPhiDirection",        "Gain on interception angle"						, &cfg_.C<St_TpcPhiDirectionC>(), 0, 0},
    {"TpcTanL",                "Gain on Tan(lambda)"							, &cfg_.C<St_TpcTanLC>(), 0, 0},
    {"TpcdXCorrectionB",       "dX correction"								, &cfg_.C<St_TpcdXCorrectionBC>(), 0, 0},
    {"TpcEffectivedX",         "dEdx correction wrt Bichsel parameterization"				, &cfg_.C<St_TpcEffectivedXC>(), 0, 0},
    {"TpcPadTBins",            "Variation on cluster size"						, 0, 0, 0},
    {"TpcZDC"        ,         "Gain on Zdc CoincidenceRate"						, &cfg_.C<St_TpcZDCC>(), 0, 0},
    {"TpcPadCorrectionMDF",    "Gain Variation along the anode wire"					, &cfg_.C<St_TpcPadCorrectionMDF>(), 0, 0},
    {"Final"        ,          ""									, 0, 0, 0},
    {"TpcLengthCorrectionB",   "Variation vs Track length and relative error in Ionization"		, &cfg_.C<St_TpcLengthCorrectionBC>(), 0, 0},
    {"TpcLengthCorrectionMDF", "Variation vs Track length and <log2(dX)> and rel. error in dE/dx"	, &cfg_.C<St_TpcLengthCorrectionMDF>(), 0, 0},
    {"TpcNoAnodeVGainC",       "Remove tpc Anode Voltage gain correction"				, 0, 0, 0},
    {"TpcdEdxCor",             "dEdx correction wrt Bichsel parameterization"                           , 0, 0, 0}
  }
{
  const St_tpcCorrectionC* chair = 0;
  const St_MDFCorrectionC* chairMDF = 0;
  const tpcCorrection* cor = 0;
  const MDFCorrection* corMDF = 0;
  int N = 0;
  int npar = 0;

  for (int k = kUncorrected + 1; k < kTpcAllCorrections; k++) {
    int nrows = 0;

    if (!corrections_[k].Chair) continue;

    if (!TESTBIT(options_, k) || corrections_[k].Chair->IsMarked()) {
      // correction is missing
      goto CLEAR;
    }

    chair    = dynamic_cast<St_tpcCorrectionC*>(corrections_[k].Chair);
    chairMDF = dynamic_cast<St_MDFCorrectionC*>(corrections_[k].Chair);

    if (! chair && ! chairMDF) {
      // This correction is not tpcCorrection or MDFCorrection type
      corrections_[k].nrows = corrections_[k].Chair->GetNRows();
      continue; // not St_tpcCorrectionC
    }

    npar = 0;

    if (chair) {

      cor = chair->Struct();
      N = chair->GetNRows();

      if (! cor || ! N) {
        goto EMPTY;
      }

      N = cor->nrows;

      for (int i = 0; i < N; i++, cor++) {
        if (cor->nrows == 0 && cor->idx == 0) continue;

        if (std::abs(cor->npar) > 0       ||
            std::abs(cor->OffSet) > 1.e-7 ||
            std::abs(cor->min)    > 1.e-7 ||
            std::abs(cor->max)    > 1.e-7) {
          npar++;
          nrows++;
        }
      }

      if (!npar) {
        // no significant corrections => switch it off
        goto CLEAR;
      }

      corrections_[k].nrows = nrows;
      continue;
    }

    corMDF = chairMDF->Struct();
    N = chairMDF->GetNRows();

    if (!corMDF || !N) {
      goto EMPTY;
    }

    npar = 0;

    for (int i = 0; i < N; i++, corMDF++) {
      if (corMDF->nrows == 0 && corMDF->idx == 0) continue;

      npar++;
      nrows++;
    }

    if (!npar) {
      // no significant corrections => switch it off
      goto CLEAR;
    }

    corrections_[k].nrows = nrows;
    continue;
EMPTY:
    // correction is empty
CLEAR:
    CLRBIT(options_, k);
    corrections_[k].Chair = 0;
  }

  // Use only one ADC correction
  if (corrections_[kAdcCorrection         ].Chair &&
      corrections_[kAdcCorrectionMDF      ].Chair) {
    LOG_ERROR << " \tTwo ADC corrections activated ? Keep active only AdcCorrectionMDF\n";
    corrections_[kAdcCorrection         ].Chair = 0;
  }
}


int StTpcdEdxCorrection::dEdxCorrection(dEdxY2_t &CdEdx)
{
  if (CdEdx.F.dE <= 0.) CdEdx.F.dE = 1;

  double dEU   = CdEdx.F.dE;
  double dE    = dEU;
  int sector   = CdEdx.sector;
  int row      = CdEdx.row;
  double dx    = CdEdx.F.dx;
  double adcCF = CdEdx.adc;

  if (dx <= 0 || (dEU <= 0 && adcCF <= 0)) return 3;

  CdEdx.channel = cfg_.C<St_TpcAvgPowerSupplyC>().ChannelFromRow(sector, row);

  CdEdx.Voltage = cfg_.C<St_tpcAnodeHVavgC>().voltagePadrow(sector, row);
  CdEdx.Crow    = cfg_.C<St_TpcAvgCurrentC>().AvCurrRow(sector, row);
  double Qcm    = cfg_.C<St_TpcAvgCurrentC>().AcChargeRowL(sector, row); // C/cm
  CdEdx.Qcm     = 1e6 * Qcm; // uC/cm

  double ZdriftDistance = CdEdx.ZdriftDistance;
  ESector kTpcOutIn = kTpcOuter;

  if (row <= cfg_.C<St_tpcPadConfigC>().innerPadRows(sector)) kTpcOutIn = kTpcInner;

  const tss_tsspar& tsspar = cfg_.S<tss_tsspar>();
  float gasGain = 1;

  if (row > cfg_.C<St_tpcPadConfigC>().innerPadRows(sector)) {
    gasGain = tpcrs::GainCorrection(sector, row, cfg_) * tsspar.wire_coupling_out;
  }
  else {
    gasGain = tpcrs::GainCorrection(sector, row, cfg_) * tsspar.wire_coupling_in;
  }

  if (gasGain <= 0.0) return 4;

  double Adc2GeVReal = tsspar.ave_ion_pot * tsspar.scale / gasGain;
  const tpcGas& tpc_gas = cfg_.S<tpcGas>();
  double ZdriftDistanceO2 = ZdriftDistance * tpc_gas.ppmOxygenIn;
  double ZdriftDistanceO2W = ZdriftDistanceO2 * tpc_gas.ppmWaterOut;
  CdEdx.ZdriftDistanceO2 = ZdriftDistanceO2;
  CdEdx.ZdriftDistanceO2W = ZdriftDistanceO2W;
  double gc, ADC, xL2, dXCorr;
  double iCut = 0;
  double slope = 0;
  int nrows = 0;
  double VarXs[kTpcLast] = {-999.};
  VarXs[kTpcZDC]               = (CdEdx.Zdc > 0) ? std::log10(CdEdx.Zdc) : 0;
  VarXs[kTpcCurrentCorrection] = CdEdx.Crow;
  VarXs[kTpcrCharge]           = CdEdx.rCharge;
  VarXs[kTpcRowQ]              = CdEdx.Qcm;
  VarXs[ktpcPressure]          = std::log(tpc_gas.barometricPressure);
  VarXs[kDrift]                = ZdriftDistanceO2;      // Blair correction
  VarXs[kMultiplicity]         = CdEdx.QRatio;
  VarXs[kzCorrection]          = ZdriftDistance;
  VarXs[ktpcMethaneIn]         = tpc_gas.percentMethaneIn * 1000. / tpc_gas.barometricPressure;
  VarXs[ktpcGasTemperature]    = tpc_gas.outputGasTemperature;
  VarXs[ktpcWaterOut]          = tpc_gas.ppmWaterOut;
  VarXs[kEdge]                 = CdEdx.PhiR;

  if (VarXs[kEdge] < -1) VarXs[kEdge] = -1;
  if (VarXs[kEdge] >  1) VarXs[kEdge] =  1;

  VarXs[kPhiDirection]     = (std::abs(CdEdx.xyzD[0]) > 1.e-7) ? std::abs(CdEdx.xyzD[1] / CdEdx.xyzD[0]) : 999.;
  VarXs[kTanL]             = CdEdx.TanL;
  VarXs[ktpcTime]          = CdEdx.tpcTime;
  VarXs[kAdcCorrection]    = adcCF;
  VarXs[kAdcCorrectionMDF] = adcCF;

  for (int k = kUncorrected; k <= kTpcLast; k++) {
    int l = 0;
    tpcCorrection* cor = 0;
    tpcCorrection* corl = 0;

    if (CdEdx.lSimulated) {
      if (k == kAdcCorrection) dE *= 2.116;//  1.75; // 1.25 is in Trs already <<<< !!!!!!

      goto ENDL;
    }

    if (! TESTBIT(options_, k)) goto ENDL;

    if (! corrections_[k].Chair) goto ENDL;

    if (k == kTpcSecRowB || k == kTpcSecRowC ) {
      const St_TpcSecRowCorC* chair = (const St_TpcSecRowCorC*) corrections_[k].Chair;

      if (! chair) goto ENDL;

      const TpcSecRowCor* gain = chair->Struct(sector - 1);
      gc =  gain->GainScale[row - 1];

      if (gc <= 0.0) return 1;

      dE *= gc;
      CdEdx.Weight = 1;

      if (gain->GainRms[row - 1] > 0.1) CdEdx.Weight = 1. / (gain->GainRms[row - 1] * gain->GainRms[row - 1]);

      goto ENDL;
    }
    else if (k == kTpcEffectivedX) {
      if      (kTpcOutIn == kTpcOuter) dx *= ((const St_TpcEffectivedXC* ) corrections_[k].Chair)->scaleOuter();
      else if (kTpcOutIn == kTpcInner ||
               kTpcOutIn == kiTpc )    dx *= ((const St_TpcEffectivedXC* ) corrections_[k].Chair)->scaleInner();

      goto ENDL;
    }

    if (k == kTpcPadMDF) {
      l = 2 * (sector - 1);

      if (row <= cfg_.C<St_tpcPadConfigC>().innerPadRows(sector)) l += kTpcInner; // for both tpc and iTPC inner sectors

      dE *= std::exp(-((St_TpcPadCorrectionMDF*)corrections_[k].Chair)->Eval(l, CdEdx.yrow, CdEdx.xpad));
      goto ENDL;
    }

    if (k == kAdcCorrectionMDF) {
      ADC = adcCF;

      if (ADC <= 0) return 3; //HACK to avoid FPE (VP)

      double xx[2] = {std::log(ADC), (double)(CdEdx.npads + CdEdx.ntmbks)};
      l = kTpcOutIn;
      int nrows = ((St_TpcAdcCorrectionMDF*) corrections_[k].Chair)->nrows();

      if (l >= nrows) l = nrows - 1;

      dE = ADC * Adc2GeVReal * ((St_TpcAdcCorrectionMDF*) corrections_[k].Chair)->Eval(l, xx);
      goto ENDL;
    }

    cor = ((St_tpcCorrectionC*) corrections_[k].Chair)->Struct();

    if (! cor) goto ENDL;

    nrows = cor->nrows;

    if (nrows <= 3)
      l = std::min(nrows - 1, static_cast<int>(kTpcOutIn));
    else {
      if (nrows == cfg_.S<tpcPadPlanes>().padRows) l = row - 1;
      else if (nrows == 192) {l = 8 * (sector - 1) + CdEdx.channel - 1; assert(l == (cor + l)->idx - 1);}
      else if (nrows ==  48) {l = 2 * (sector - 1) + kTpcOutIn;}
      else if (nrows ==   6) {l =                    kTpcOutIn;             if (sector > 12) l += 3;}
      else if (nrows ==   4) {l = std::min(static_cast<int>(kTpcOutIn), 1); if (sector > 12) l += 2;}
    }

    corl = cor + l;
    iCut = 0;

    if (k == kAdcCorrection) {
      ADC = adcCF;

      if (ADC <= 0) return 3; //HACK to avoid FPE (VP)

      if (corl->type == 12)
        dE = Adc2GeVReal * ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(l, ADC, VarXs[kTanL]);
      else
        dE = Adc2GeVReal * ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(l, ADC, std::abs(CdEdx.zG));

      if (dE <= 0) return 3;

      goto ENDL;
    }
    else if (k == kTpcdCharge) {
      if (l > 2) l = 1;

      slope = ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(l, row + 0.5);
      dE *=  std::exp(-slope * CdEdx.dCharge);
      dE *=  std::exp(-((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(2 + l, CdEdx.dCharge));
      goto ENDL;
    }
    else if (k == kzCorrection) {
      iCut = 1; // Always cut
    }
    else if (k == kdXCorrection) {
      xL2 = std::log2(dx);
      dXCorr = ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(l, xL2);

      if (std::abs(dXCorr) > 10) return 3;

      if (nrows == 7) {// old schema without iTPC
        dXCorr += ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(2, xL2);
        dXCorr += ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(5 + kTpcOutIn, xL2);
      }

      CdEdx.dxC = std::exp(dXCorr) * CdEdx.F.dx;
      dE *= std::exp(-dXCorr);
      goto ENDL;
    }
    else if (k == kSpaceCharge) {
      if (cor[2 * l  ].min <= CdEdx.QRatio && CdEdx.QRatio <= cor[2 * l  ].max &&
          cor[2 * l + 1].min <= CdEdx.DeltaZ && CdEdx.DeltaZ <= cor[2 * l + 1].max)
        dE *= std::exp(-((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(2 * l, CdEdx.QRatio)
                         - ((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(2 * l + 1, CdEdx.DeltaZ));

      goto ENDL;
    }
    else if (k == kEdge) {
      if (corl->type == 200) VarXs[kEdge] = std::abs(CdEdx.edge);

      if (corl->min > 0 && corl->min > VarXs[kEdge]    ) return 2;
    }
    else if (k == ktpcTime) {   // use the correction if you have xmin < xmax && xmin <= x <= xmax
      if (corl->min >= corl->max || corl->min > VarXs[ktpcTime] ||  VarXs[ktpcTime] > corl->max) goto ENDL;

      double xx = VarXs[ktpcTime];
      dE *= std::exp(-((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(l, xx));
      goto ENDL;
    }

    if (corl->type == 300) {
      if (corl->min > 0 && corl->min > VarXs[k]    ) VarXs[k] = corl->min;

      if (corl->max > 0 && VarXs[k]     > corl->max) VarXs[k] = corl->max;
    }

    if (std::abs(corl->npar) >= 100 || iCut) {
      int iok = 2;

      if (corl->min >= corl->max) {
        iok = 0;
      }
      else {
        if (corl->min <= VarXs[k] && VarXs[k] <= corl->max) {
          iok = 0;
        }
      }

      if (iok) {
        return iok;
      }
    }

    if (corl->npar % 100) {
      double dECor = std::exp(-((St_tpcCorrectionC*)corrections_[k].Chair)->CalcCorrection(l, VarXs[k]));
      dE *= dECor;
    }

ENDL:
    CdEdx.C[k].dE = dE;
    CdEdx.C[k].dx = dx;
    CdEdx.C[k].dEdx  = dE/dx;
    CdEdx.C[k].dEdxL = std::log(CdEdx.C[k].dEdx);
  }

  CdEdx.F = CdEdx.C[kTpcLast];
  return 0;
}


void StTpcdEdxCorrection::Print(const dEdxY2_t& mdEdx) const
{
  LOG_INFO << "StTpcdEdxCorrection:: Sector/row/pad " << mdEdx.sector << "/" << mdEdx.row << "/" << mdEdx.pad << '\n'
           << "\tdrift distance / O2 / O2W " << mdEdx.ZdriftDistance << "/" << mdEdx.ZdriftDistanceO2 << "/" << mdEdx.ZdriftDistanceO2W << '\n'
           << "Local xyz " << mdEdx.xyz[0] << "\t" << mdEdx.xyz[1] << "\t" << mdEdx.xyz[2] << '\n'
           << "Local xyzD " << mdEdx.xyzD[0] << "\t" << mdEdx.xyzD[1] << "\t" << mdEdx.xyzD[2] << '\n';

  for (int k = (int)kUncorrected; k <= ((int)kTpcLast) + 1; k++)
  {
    LOG_INFO << k;

    if (k <= (int) kTpcLast) {
      LOG_INFO << "\tdE "         << mdEdx.C[k].dE
               << "\tdx "         << mdEdx.C[k].dx
               << "\tdE/dx "      << mdEdx.C[k].dEdx
               << "\tlog(dE/dx) " << mdEdx.C[k].dEdxL
               << '\t'            << corrections_[k].name
               << '\t'            << corrections_[k].title;
    }
    else {
      LOG_INFO << "\tdE "         << mdEdx.F.dE
               << "\tdx "         << mdEdx.F.dx
               << "\tdE/dx "      << mdEdx.F.dEdx
               << "\tlog(dE/dx) " << mdEdx.F.dEdxL
               << "\tFinal \t ";
    }

    LOG_INFO << '\n';
  }
}
