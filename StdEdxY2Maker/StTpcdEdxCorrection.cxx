/*!
  \class dEdxY2_t
  dEdxY2_t class contains data for each cluster which is used during different calibartion steps
  \class StTpcdEdxCorrection
  StTpcdEdxCorrection class contains utilities
*/

#include <iostream>

#include "StdEdxY2Maker/StTpcdEdxCorrection.h"
#include "StTpcDb/StTpcDb.h"
#include "StBichsel/Bichsel.h"
#include "StDetectorDbMaker/St_tss_tssparC.h"
#include "StDetectorDbMaker/St_TpcEdgeC.h"
#include "StDetectorDbMaker/St_TpcAdcCorrectionBC.h"
#include "StDetectorDbMaker/St_TpcAdcCorrectionMDF.h"
#include "StDetectorDbMaker/St_TpcdChargeC.h"
#include "StDetectorDbMaker/St_TpcrChargeC.h"
#include "StDetectorDbMaker/St_TpcCurrentCorrectionC.h"
#include "StDetectorDbMaker/St_TpcRowQC.h"
#include "StDetectorDbMaker/St_TpcSecRowBC.h"
#include "StDetectorDbMaker/St_TpcSecRowCC.h"
#include "StDetectorDbMaker/St_tpcPressureBC.h"
#include "StDetectorDbMaker/St_TpcDriftDistOxygenC.h"
#include "StDetectorDbMaker/St_TpcMultiplicityC.h"
#include "StDetectorDbMaker/St_TpcZCorrectionBC.h"
#include "StDetectorDbMaker/St_tpcMethaneInC.h"
#include "StDetectorDbMaker/St_tpcGasC.h"
#include "StDetectorDbMaker/St_tpcGasTemperatureC.h"
#include "StDetectorDbMaker/St_tpcWaterOutC.h"
#include "StDetectorDbMaker/St_TpcSpaceChargeC.h"
#include "StDetectorDbMaker/St_TpcPhiDirectionC.h"
#include "StDetectorDbMaker/St_TpcTanLC.h"
#include "StDetectorDbMaker/St_TpcdXCorrectionBC.h"
#include "StDetectorDbMaker/St_TpcEffectivedXC.h"
#include "StDetectorDbMaker/St_TpcZDCC.h"
#include "StDetectorDbMaker/St_TpcLengthCorrectionBC.h"
#include "StDetectorDbMaker/St_TpcLengthCorrectionMDF.h"
#include "StDetectorDbMaker/St_TpcPadCorrectionMDF.h"
#include "StDetectorDbMaker/St_TpcdEdxCorC.h"
#include "StDetectorDbMaker/St_tpcAnodeHVavgC.h"
#include "StDetectorDbMaker/St_TpcAvgCurrentC.h"
#include "StDetectorDbMaker/St_TpcAvgPowerSupplyC.h"
#include "StDetectorDbMaker/St_tpcTimeDependenceC.h"
#include "StDetectorDbMaker/St_trigDetSumsC.h"
#include "tpcrs/logger.h"


StTpcdEdxCorrection::StTpcdEdxCorrection(int option, int debug) :
  m_Mask(option), m_tpcGas(0),
  m_Debug(debug)
{
  assert(gStTpcDb);

  if (!m_Mask) m_Mask = -1;

  static const char* FXTtables[] = {"TpcdXCorrectionB",
                                      "tpcGainCorrection",
                                      "TpcLengthCorrectionMDF",
                                      "TpcPadCorrectionMDF",
                                      "TpcSecRowB",
                                      "TpcZCorrectionB"
                                     };
  static int NT = sizeof(FXTtables) / sizeof(const char*);

  ReSetCorrections();
}


void StTpcdEdxCorrection::ReSetCorrections()
{
  St_tpcGasC* tpcGas = (St_tpcGasC*) St_tpcGasC::instance();  //

  if (!tpcGas || ! tpcGas->GetNRows()) {
    LOG_ERROR << "=== tpcGas is missing ===\n";
    assert(tpcGas);
  }

  SettpcGas(tpcGas);
  memset (m_Corrections, 0, sizeof(m_Corrections));
  m_Corrections[kUncorrected           ] = dEdxCorrection_t("UnCorrected",            ""								, 0);
  m_Corrections[kAdcCorrection         ] = dEdxCorrection_t("TpcAdcCorrectionB",      "ADC/Clustering nonlinearity correction"				, St_TpcAdcCorrectionBC::instance());
  m_Corrections[kEdge                  ] = dEdxCorrection_t("TpcEdge",                "Gain on distance from Chamber edge"				, St_TpcEdgeC::instance());
  m_Corrections[kAdcCorrectionMDF      ] = dEdxCorrection_t("TpcAdcCorrectionMDF",    "ADC/Clustering nonlinearity correction MDF"			, St_TpcAdcCorrectionMDF::instance());
  m_Corrections[kTpcdCharge            ] = dEdxCorrection_t("TpcdCharge",             "ADC/Clustering undershoot correction"				, St_TpcdChargeC::instance());
  m_Corrections[kTpcrCharge            ] = dEdxCorrection_t("TpcrCharge",             "ADC/Clustering rounding correction"				, St_TpcrChargeC::instance());
  m_Corrections[kTpcCurrentCorrection  ] = dEdxCorrection_t("TpcCurrentCorrection",   "Correction due to sagg of Voltage due to anode current"		, St_TpcCurrentCorrectionC::instance());
  m_Corrections[kTpcRowQ               ] = dEdxCorrection_t("TpcRowQ",                "Gas gain correction for row versus accumulated charge,"		, St_TpcRowQC::instance());
  m_Corrections[kTpcSecRowB            ] = dEdxCorrection_t("TpcSecRowB",             "Gas gain correction for sector/row"				, St_TpcSecRowBC::instance());
  m_Corrections[kTpcSecRowC            ] = dEdxCorrection_t("TpcSecRowC",             "Additional Gas gain correction for sector/row"			, St_TpcSecRowCC::instance());
  m_Corrections[ktpcPressure           ] = dEdxCorrection_t("tpcPressureB",           "Gain on Gas Density due to Pressure"				, St_tpcPressureBC::instance());
  m_Corrections[ktpcTime               ] = dEdxCorrection_t("tpcTime"       ,         "Unregognized time dependce"					, St_tpcTimeDependenceC::instance());
  m_Corrections[kDrift                 ] = dEdxCorrection_t("TpcDriftDistOxygen",     "Correction for Electron Attachment due to O2"			, St_TpcDriftDistOxygenC::instance());
  m_Corrections[kMultiplicity          ] = dEdxCorrection_t("TpcMultiplicity",        "Global track multiplicity dependence"				, St_TpcMultiplicityC::instance());
  m_Corrections[kzCorrection           ] = dEdxCorrection_t("TpcZCorrectionB",        "Variation on drift distance"					, St_TpcZCorrectionBC::instance());
  m_Corrections[ktpcMethaneIn          ] = dEdxCorrection_t("tpcMethaneIn",           "Gain on Methane content"						, St_tpcMethaneInC::instance());
  m_Corrections[ktpcGasTemperature     ] = dEdxCorrection_t("tpcGasTemperature",      "Gain on Gas Dens. due to Temperature"				, St_tpcGasTemperatureC::instance());
  m_Corrections[ktpcWaterOut           ] = dEdxCorrection_t("tpcWaterOut",            "Gain on Water content"						, St_tpcWaterOutC::instance());
  m_Corrections[kSpaceCharge           ] = dEdxCorrection_t("TpcSpaceCharge",         "Gain on space charge near the wire"				, St_TpcSpaceChargeC::instance());
  m_Corrections[kPhiDirection          ] = dEdxCorrection_t("TpcPhiDirection",        "Gain on interception angle"					, St_TpcPhiDirectionC::instance());
  m_Corrections[kTanL                  ] = dEdxCorrection_t("TpcTanL",                "Gain on Tan(lambda)"						, St_TpcTanLC::instance());
  m_Corrections[kdXCorrection          ] = dEdxCorrection_t("TpcdXCorrectionB",       "dX correction"							, St_TpcdXCorrectionBC::instance());
  m_Corrections[kTpcEffectivedX        ] = dEdxCorrection_t("TpcEffectivedX",         "dEdx correction wrt Bichsel parameterization"			, St_TpcEffectivedXC::instance());
  m_Corrections[kTpcPadTBins           ] = dEdxCorrection_t("TpcPadTBins",            "Variation on cluster size"					, 0);
  m_Corrections[kTpcZDC                ] = dEdxCorrection_t("TpcZDC"        ,         "Gain on Zdc CoincidenceRate"					, St_TpcZDCC::instance());
  m_Corrections[kTpcPadMDF             ] = dEdxCorrection_t("TpcPadCorrectionMDF",    "Gain Variation along the anode wire"				, St_TpcPadCorrectionMDF::instance());
  m_Corrections[kTpcLast               ] = dEdxCorrection_t("Final"        ,          ""								, 0);
  m_Corrections[kTpcLengthCorrection   ] = dEdxCorrection_t("TpcLengthCorrectionB",   "Variation vs Track length and relative error in Ionization"	, St_TpcLengthCorrectionBC::instance());
  m_Corrections[kTpcLengthCorrectionMDF] = dEdxCorrection_t("TpcLengthCorrectionMDF", "Variation vs Track length and <log2(dX)> and rel. error in dE/dx", St_TpcLengthCorrectionMDF::instance());
  m_Corrections[kTpcNoAnodeVGainC      ] = dEdxCorrection_t("TpcNoAnodeVGainC",       "Remove tpc Anode Voltage gain correction"			, 0);
  m_Corrections[kTpcdEdxCor            ] = dEdxCorrection_t("TpcdEdxCor",             "dEdx correction wrt Bichsel parameterization"			, St_TpcdEdxCorC::instance());
  const St_tpcCorrectionC* chair = 0;
  const St_MDFCorrectionC* chairMDF = 0;
  const tpcCorrection_st* cor = 0;
  const MDFCorrection_st* corMDF = 0;
  int N = 0;
  int npar = 0;
  int nrows = 0;

  for (int k = kUncorrected + 1; k < kTpcAllCorrections; k++) {
    if (! m_Corrections[k].Chair) continue;

    nrows = 0;
    LOG_INFO << "StTpcdEdxCorrection: " << m_Corrections[k].Name << "/" << m_Corrections[k].Title << '\n';

    if (! TESTBIT(m_Mask, k) || m_Corrections[k].Chair->IsMarked()) {
      LOG_INFO << " \tis missing\n";
      goto CLEAR;
    }

    chair    = dynamic_cast<St_tpcCorrectionC*>(m_Corrections[k].Chair);
    chairMDF = dynamic_cast<St_MDFCorrectionC*>(m_Corrections[k].Chair);

    if (! chair && ! chairMDF) {
      LOG_WARN << " \tis not tpcCorrection or MDFCorrection type\n";
      m_Corrections[k].nrows = m_Corrections[k].Chair->GetNRows();
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

      if (! npar ) {
        LOG_INFO << " \thas no significant corrections => switch it off\n";
        goto CLEAR;
      }

      m_Corrections[k].nrows = nrows;
      continue;
    }


    corMDF = chairMDF->Struct();
    N = chairMDF->GetNRows();

    if (! corMDF || ! N) {
      goto EMPTY;
    }

    npar = 0;

    for (int i = 0; i < N; i++, corMDF++) {
      if (corMDF->nrows == 0 && corMDF->idx == 0) continue;

      npar++;
      nrows++;
    }

    if (! npar ) {
      LOG_INFO << " \thas no significant corrections => switch it off\n";
      goto CLEAR;
    }

    m_Corrections[k].nrows = nrows;
    continue;
EMPTY:
    LOG_INFO << " \tis empty\n";
CLEAR:
    CLRBIT(m_Mask, k);
    m_Corrections[k].Chair = 0;
  }

  // Use only one ADC correction
  if (m_Corrections[kAdcCorrection         ].Chair &&
      m_Corrections[kAdcCorrectionMDF      ].Chair) {
    LOG_ERROR << " \tTwo ADC corrections activated ? Keep active only AdcCorrectionMDF\n";
    m_Corrections[kAdcCorrection         ].Chair = 0;
  }
}


StTpcdEdxCorrection::~StTpcdEdxCorrection()
{
  // Can't delete because the chairs are also used in StTpcRSMaker
  //  for (int k = 0; k < kTpcAllCorrections; k++) SafeDelete(m_Corrections[k].Chair);
}


int  StTpcdEdxCorrection::dEdxCorrection(dEdxY2_t &CdEdx, bool doIT)
{
  //  static const double Degree2Rad = TMath::Pi()/180.;
  mdEdx = &CdEdx;

  if (CdEdx.F.dE <= 0.) CdEdx.F.dE = 1;

  double dEU = CdEdx.F.dE;
  double dE  = dEU;
  int sector            = CdEdx.sector;
  int row       	  = CdEdx.row;
  double dx     	  = CdEdx.F.dx;
  double adcCF = CdEdx.adc;

  if (dx <= 0 || (dEU <= 0 && adcCF <= 0)) return 3;

  int channel = St_TpcAvgPowerSupplyC::instance()->ChannelFromRow(sector, row);
  CdEdx.channel = channel;

  CdEdx.Voltage = St_tpcAnodeHVavgC::instance()->voltagePadrow(sector, row);
  CdEdx.Crow    = St_TpcAvgCurrentC::instance()->AvCurrRow(sector, row);
  double    Qcm      = St_TpcAvgCurrentC::instance()->AcChargeRowL(sector, row); // C/cm
  CdEdx.Qcm     = 1e6 * Qcm; // uC/cm

  double ZdriftDistance = CdEdx.ZdriftDistance;
  ESector kTpcOutIn = kTpcOuter;

  if (! St_tpcPadConfigC::instance()->iTpc(sector)) {
    if (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) kTpcOutIn = kTpcInner;
  }
  else {
    if (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) kTpcOutIn = kiTpc;
  }

  St_tss_tssparC* tsspar = St_tss_tssparC::instance();
  float gasGain = 1;
  float gainNominal = 0;

  if (row > St_tpcPadConfigC::instance()->innerPadRows(sector)) {
    gainNominal = tsspar->gain_out() * tsspar->wire_coupling_out();
    gasGain = tsspar->gain_out(sector, row) * tsspar->wire_coupling_out();
  }
  else {
    gainNominal = tsspar->gain_in() * tsspar->wire_coupling_in();
    gasGain = tsspar->gain_in(sector, row) * tsspar->wire_coupling_in();
  }

  if (gasGain <= 0.0) return 4;

  //  double gainAVcorr = gasGain/gainNominal;
  mAdc2GeV = tsspar->ave_ion_pot() * tsspar->scale() / gainNominal;
  double Adc2GeVReal = tsspar->ave_ion_pot() * tsspar->scale() / gasGain;
  tpcGas_st* gas = m_tpcGas->Struct();
  double ZdriftDistanceO2 = ZdriftDistance * gas->ppmOxygenIn;
  double ZdriftDistanceO2W = ZdriftDistanceO2 * gas->ppmWaterOut;
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
  VarXs[kTpcPadTBins]          = CdEdx.Npads * CdEdx.Ntbins;
  VarXs[ktpcPressure]          = std::log(gas->barometricPressure);
  VarXs[kDrift]                = ZdriftDistanceO2;      // Blair correction
  VarXs[kMultiplicity]         = CdEdx.QRatio;
  VarXs[kzCorrection]          = ZdriftDistance;
  VarXs[ktpcMethaneIn]         = gas->percentMethaneIn * 1000. / gas->barometricPressure;
  VarXs[ktpcGasTemperature]    = gas->outputGasTemperature;
  VarXs[ktpcWaterOut]          = gas->ppmWaterOut;
  VarXs[kEdge]                 = CdEdx.PhiR;

  if (VarXs[kEdge] < -1) VarXs[kEdge] = -1;

  if (VarXs[kEdge] >  1) VarXs[kEdge] =  1;

  VarXs[kPhiDirection]         = (std::abs(CdEdx.xyzD[0]) > 1.e-7) ? std::abs(CdEdx.xyzD[1] / CdEdx.xyzD[0]) : 999.;
  VarXs[kTanL]                 = CdEdx.TanL;
  VarXs[ktpcTime]              = CdEdx.tpcTime;
  VarXs[kAdcCorrection] = VarXs[kAdcCorrectionMDF] = adcCF;

  for (int k = kUncorrected; k <= kTpcLast; k++) {
    int l = 0;
    tpcCorrection_st* cor = 0;
    tpcCorrection_st* corl = 0;

    if (CdEdx.lSimulated) {
      if (k == kAdcCorrection) dE *= 2.116;//  1.75; // 1.25 is in Trs already <<<< !!!!!!

      goto ENDL;
    }

    if (! TESTBIT(m_Mask, k)) goto ENDL;

    if (! m_Corrections[k].Chair) goto ENDL;

    if (k == kTpcSecRowB || k == kTpcSecRowC ) {
      const St_TpcSecRowCorC* chair = (const St_TpcSecRowCorC*) m_Corrections[k].Chair;

      if (! chair) goto ENDL;

      const TpcSecRowCor_st* gain = chair->Struct(sector - 1);
      gc =  gain->GainScale[row - 1];

      if (gc <= 0.0) return 1;

      dE *= gc;
      CdEdx.Weight = 1;

      if (gain->GainRms[row - 1] > 0.1) CdEdx.Weight = 1. / (gain->GainRms[row - 1] * gain->GainRms[row - 1]);

      goto ENDL;
    }
    else if (k == kTpcEffectivedX) {
      if      (kTpcOutIn == kTpcOuter) dx *= ((const St_TpcEffectivedXC* ) m_Corrections[k].Chair)->scaleOuter();
      else if (kTpcOutIn == kTpcInner ||
               kTpcOutIn == kiTpc )    dx *= ((const St_TpcEffectivedXC* ) m_Corrections[k].Chair)->scaleInner();

      goto ENDL;
    }

    if (k == kTpcPadMDF) {
      l = 2 * (sector - 1);

      if (row <= St_tpcPadConfigC::instance()->innerPadRows(sector)) l += kTpcInner; // for both tpc and iTPC inner sectors

      dE *= std::exp(-((St_TpcPadCorrectionMDF*)m_Corrections[k].Chair)->Eval(l, CdEdx.yrow, CdEdx.xpad));
      goto ENDL;
    }

    if (k == kAdcCorrectionMDF) {
      ADC = adcCF;

      if (ADC <= 0) return 3; //HACK to avoid FPE (VP)

      double xx[2] = {std::log(ADC), (double)(CdEdx.npads + CdEdx.ntmbks)};
      l = kTpcOutIn;
      int nrows = ((St_TpcAdcCorrectionMDF*) m_Corrections[k].Chair)->nrows();

      if (l >= nrows) l = nrows - 1;

      dE = ADC * Adc2GeVReal * ((St_TpcAdcCorrectionMDF*) m_Corrections[k].Chair)->Eval(l, xx);
      goto ENDL;
    }

    cor = ((St_tpcCorrectionC*) m_Corrections[k].Chair)->Struct();

    if (! cor) goto ENDL;

    nrows = cor->nrows;

    if (nrows <= 3)
      l = std::min(nrows - 1, static_cast<int>(kTpcOutIn));
    else {
      if (nrows == St_tpcPadConfigC::instance()->numberOfRows(sector)) l = row - 1;
      else if (nrows == 192) {l = 8 * (sector - 1) + channel - 1; assert(l == (cor + l)->idx - 1);}
      else if (nrows ==  48) {l = 2 * (sector - 1) + kTpcOutIn;}
      else if (nrows ==   6) {l =            kTpcOutIn;     if (sector > 12) l += 3;}
      else if (nrows ==   4) {l = std::min(static_cast<int>(kTpcOutIn), 1); if (sector > 12) l += 2;}
    }

    corl = cor + l;
    iCut = 0;

    if (k == kAdcCorrection) {
      ADC = adcCF;

      if (ADC <= 0) return 3; //HACK to avoid FPE (VP)

      if (corl->type == 12)
        dE = Adc2GeVReal * ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(l, ADC, VarXs[kTanL]);
      else
        dE = Adc2GeVReal * ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(l, ADC, std::abs(CdEdx.zG));

      if (dE <= 0) return 3;

      goto ENDL;
    }
    else if (k == kTpcdCharge) {
      if (l > 2) l = 1;

      slope = ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(l, row + 0.5);
      dE *=  std::exp(-slope * CdEdx.dCharge);
      dE *=  std::exp(-((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(2 + l, CdEdx.dCharge));
      goto ENDL;
    }
    else if (k == kzCorrection) {
      iCut = 1; // Always cut
    }
    else if (k == kdXCorrection) {
      xL2 = std::log2(dx);
      dXCorr = ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(l, xL2);

      if (std::abs(dXCorr) > 10) return 3;

      if (nrows == 7) {// old schema without iTPC
        dXCorr += ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(2, xL2);
        dXCorr += ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(5 + kTpcOutIn, xL2);
      }

      CdEdx.dxC = std::exp(dXCorr) * CdEdx.F.dx;
      dE *= std::exp(-dXCorr);
      goto ENDL;
    }
    else if (k == kSpaceCharge) {
      if (cor[2 * l  ].min <= CdEdx.QRatio && CdEdx.QRatio <= cor[2 * l  ].max &&
          cor[2 * l + 1].min <= CdEdx.DeltaZ && CdEdx.DeltaZ <= cor[2 * l + 1].max)
        dE *= std::exp(-((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(2 * l, CdEdx.QRatio)
                         - ((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(2 * l + 1, CdEdx.DeltaZ));

      goto ENDL;
    }
    else if (k == kEdge) {
      if (corl->type == 200) VarXs[kEdge] = std::abs(CdEdx.edge);

      if (corl->min > 0 && corl->min > VarXs[kEdge]    ) return 2;
    }
    else if (k == ktpcTime) {   // use the correction if you have xmin < xmax && xmin <= x <= xmax
      if (corl->min >= corl->max || corl->min > VarXs[ktpcTime] ||  VarXs[ktpcTime] > corl->max) goto ENDL;

      double xx = VarXs[ktpcTime];
      dE *= std::exp(-((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(l, xx));
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
      double dECor = std::exp(-((St_tpcCorrectionC*)m_Corrections[k].Chair)->CalcCorrection(l, VarXs[k]));
      dE *= dECor;
    }

ENDL:
    CdEdx.C[k].dE = dE;
    CdEdx.C[k].dx = dx;
    CdEdx.C[k].dEdx    = CdEdx.C[k].dE / CdEdx.C[k].dx;
    CdEdx.C[k].dEdxL   = std::log(CdEdx.C[k].dEdx);
  }

  CdEdx.F = CdEdx.C[kTpcLast];
  return 0;
}


void StTpcdEdxCorrection::Print(Option_t* opt) const
{
  if (! mdEdx) return;

  LOG_INFO << "StTpcdEdxCorrection:: Sector/row/pad " << mdEdx->sector << "/" << mdEdx->row << "/" << mdEdx->pad << '\n';
  LOG_INFO << "Npads/Ntbins " << mdEdx->Npads << "/" << mdEdx->Ntbins
       << "\tdrift distance / O2 / O2W " << mdEdx->ZdriftDistance << "/" << mdEdx->ZdriftDistanceO2 << "/" << mdEdx->ZdriftDistanceO2W << '\n';
  LOG_INFO << "Local xyz " << mdEdx->xyz[0] << "\t" << mdEdx->xyz[1] << "\t" << mdEdx->xyz[2] << '\n';
  LOG_INFO << "Local xyzD " << mdEdx->xyzD[0] << "\t" << mdEdx->xyzD[1] << "\t" << mdEdx->xyzD[2] << '\n';
  TString Line;

  for (int k = (int)kUncorrected; k <= ((int)kTpcLast) + 1; k++) {
    Line  = Form("%2i", k);

    if (k <= (int) kTpcLast) {
      Line += Form("\tdE %10.5g", mdEdx->C[k].dE);
      Line += Form("\tdx  %10.5g", mdEdx->C[k].dx);
      Line += Form("\tdE/dx  %10.5g", mdEdx->C[k].dEdx);
      Line += Form("\tlog(dE/dx)  %10.5g", mdEdx->C[k].dEdxL);
      Line += "\t"; Line += TString(m_Corrections[k].Name); Line += "\t"; Line +=  TString(m_Corrections[k].Title);
    }
    else {
      Line += Form("\tdE %10.5g", mdEdx->F.dE);
      Line += Form("\tdx  %10.5g", mdEdx->F.dx);
      Line += Form("\tdE/dx  %10.5g", mdEdx->F.dEdx);
      Line += Form("\tlog(dE/dx)  %10.5g", mdEdx->F.dEdxL);
      Line +=  "\tFinal \t ";
    }

    LOG_INFO << Line.Data() << '\n';
  }
}
