/* TpcResponseSimulator.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    TpcResponseSimulator.idl

  Table: TpcResponseSimulator

       description:   Tpc Response Simulator parameters  Tpc Response Simulator parameters


  2.145, effective reduction of OmegaTau near
				  Inner sector anode wire
  1.8, effective reduction of OmegaTau near
				  Outer sector anode wire
 */
#ifndef TPCRESPONSESIMULATOR_H
#define TPCRESPONSESIMULATOR_H
#define TPCRESPONSESIMULATOR_SPEC \
"struct TpcResponseSimulator { \
	float I0; \
	float Cluster; \
	float W; \
	float OmegaTau; \
	float K3IP; \
	float K3IR; \
	float K3OP; \
	float K3OR; \
	float FanoFactor; \
	float AveragePedestal; \
	float AveragePedestalRMS; \
	float AveragePedestalRMSX; \
	float tauIntegration; \
	float tauF; \
	float tauP; \
	float tauXI; \
	float tauXO; \
	float tauCI; \
	float tauCO; \
	float SigmaJitterTI; \
	float SigmaJitterTO; \
	float SigmaJitterXI; \
	float SigmaJitterXO; \
	float longitudinalDiffusion; \
	float transverseDiffusion; \
	float NoElPerAdc; \
	float NoElPerAdcI; \
	float NoElPerAdcO; \
	float NoElPerAdcX; \
	float OmegaTauScaleI; \
	float OmegaTauScaleO; \
	float SecRowCorIW[2]; \
	float SecRowCorOW[2]; \
	float SecRowCorIE[2]; \
	float SecRowCorOE[2]; \
	float SecRowSigIW[2]; \
	float SecRowSigOW[2]; \
	float SecRowSigIE[2]; \
	float SecRowSigOE[2]; \
	float PolyaInner; \
	float PolyaOuter; \
	float T0offset; \
	float T0offsetI; \
	float T0offsetO; \
	float FirstRowC; \
};"
typedef struct TpcResponseSimulator_st {
  float I0; /* = 13.1 eV, CH4 */
  float Cluster; /* = 3.2, average no. of electrons per primary  */
  float W; /* = 26.2 eV */
  float OmegaTau; /* = 3.02, fit of data */
  float K3IP; /* = 0.68,(pads) for a/s = 2.5e-3 and h/s = 0.5 */
  float K3IR; /* = 0.89,(row)  for a/s = 2.5e-3 and h/s = 0.5 */
  float K3OP; /* = 0.55,(pads) for a/s = 2.5e-3 and h/s = 1.0 */
  float K3OR; /* = 0.61,(row)  for a/s = 2.5e-3 and h/s = 1.0 */
  float FanoFactor; /* = 0.3 */
  float AveragePedestal; /* = 50.0 */
  float AveragePedestalRMS; /* = 1.4, Old Tpc electronics */
  float AveragePedestalRMSX; /* = 0.7, New Tpx electronics */
  float tauIntegration; /* = 2.5*74.6e-9  secs */
  float tauF; /* = 394.0e-9 secs Tpc */
  float tauP; /* = 775.0e-9 secs Tpc */
  float tauXI; /* =  60.0e-9 secs Tpx Inner integration time */
  float tauXO; /* =  74.6e-9  secs Tpx Outer integration time */
  float tauCI; /* =   0  */
  float tauCO; /* =   0  */
  float SigmaJitterTI; /* = 0.2  for Tpx inner */
  float SigmaJitterTO; /* = 0.2  for Tpx outer */
  float SigmaJitterXI; /* = 0.0  for Tpx inner */
  float SigmaJitterXO; /* = 0.0  for Tpx outer */
  float longitudinalDiffusion; /*   cm/sqrt(cm)  */
  float transverseDiffusion; /*   cm/sqrt(cm)  */
  float NoElPerAdc; /* = 335, No. of electrons per 1 ADC count, keep for back compartibility */
  float NoElPerAdcI; /* = 335, No. of electrons per 1 ADC count for inner TPX */
  float NoElPerAdcO; /* = 335, No. of electrons per 1 ADC count for outer TPX */
  float NoElPerAdcX; /* = 335, No. of electrons per 1 ADC count for iTPC      */
  float OmegaTauScaleI;
  float OmegaTauScaleO;
  float SecRowCorIW[2]; /* parameterization of Inner West correction vs row */
  float SecRowCorOW[2]; /* parameterization of Outer West correction vs row */
  float SecRowCorIE[2]; /* parameterization of Inner East correction vs row */
  float SecRowCorOE[2]; /* parameterization of Outer East correction vs row */
  float SecRowSigIW[2]; /* parameterization of Inner West gain sigma vs row */
  float SecRowSigOW[2]; /* parameterization of Outer West gain sigma vs row */
  float SecRowSigIE[2]; /* parameterization of Inner East gain sigma vs row */
  float SecRowSigOE[2]; /* parameterization of Outer East gain sigma vs row */
  float PolyaInner; /* = 1.38, Polya parameter for inner sectors */
  float PolyaOuter; /* = 1.38, Polya parameter for outer sectors */
  float T0offset; /* = 0.0   extra off set for Altro chip */
  float T0offsetI; /* = 0.0   extra off set for inner sector */
  float T0offsetO; /* = 0.0   extra off set for outer sector */
  float FirstRowC; /* = 0.0   extra correction for the first pad row */
} TPCRESPONSESIMULATOR_ST;
#endif /* TPCRESPONSESIMULATOR_H */
