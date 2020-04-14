#ifndef St_TpcResponseSimulatorC_h
#define St_TpcResponseSimulatorC_h

#include "tpcrs/config_structs.h"
#include "TpcResponseSimulator.h"

struct St_TpcResponseSimulatorC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcResponseSimulatorC, TpcResponseSimulator_st>
{
  float 	I0(int i = 0) 	const {return Struct(i)->I0;}
  float 	Cluster(int i = 0) 	const {return Struct(i)->Cluster;}
  float 	W(int i = 0) 	const {return Struct(i)->W;}
  float 	OmegaTau(int i = 0) 	const {return Struct(i)->OmegaTau;}
  float 	K3IP(int i = 0) 	const {return Struct(i)->K3IP;}
  float 	K3IR(int i = 0) 	const {return Struct(i)->K3IR;}
  float 	K3OP(int i = 0) 	const {return Struct(i)->K3OP;}
  float 	K3OR(int i = 0) 	const {return Struct(i)->K3OR;}
  float 	FanoFactor(int i = 0) 	const {return Struct(i)->FanoFactor;}
  float 	AveragePedestal(int i = 0) 	const {return Struct(i)->AveragePedestal;}
  float 	AveragePedestalRMS(int i = 0) 	const {return Struct(i)->AveragePedestalRMS;}
  float 	AveragePedestalRMSX(int i = 0) 	const {return Struct(i)->AveragePedestalRMSX;}
  float 	tauIntegration(int i = 0) 	const {return Struct(i)->tauIntegration;}
  float 	tauF(int i = 0) 	const {return Struct(i)->tauF;}
  float 	tauP(int i = 0) 	const {return Struct(i)->tauP;}
  float*      tauX(int i = 0) 	const {return &Struct(i)->tauXI;}
  float 	tauXI(int i = 0) 	const {return Struct(i)->tauXI;}
  float 	tauXO(int i = 0) 	const {return Struct(i)->tauXO;}
  float*      tauC(int i = 0) 	const {return &Struct(i)->tauCI;}
  float 	tauCI(int i = 0) 	const {return Struct(i)->tauCI;}
  float 	tauCO(int i = 0) 	const {return Struct(i)->tauCO;}
  float 	SigmaJitterTI(int i = 0) 	const {return Struct(i)->SigmaJitterTI;}
  float 	SigmaJitterTO(int i = 0) 	const {return Struct(i)->SigmaJitterTO;}
  float 	SigmaJitterXI(int i = 0) 	const {return Struct(i)->SigmaJitterXI;}
  float 	SigmaJitterXO(int i = 0) 	const {return Struct(i)->SigmaJitterXO;}
  float 	longitudinalDiffusion(int i = 0) 	const {return Struct(i)->longitudinalDiffusion;}
  float 	transverseDiffusion(int i = 0) 	const {return Struct(i)->transverseDiffusion;}
  float       NoElPerAdc(int i = 0)  const {return Struct(i)->NoElPerAdc;}
  float       NoElPerAdcI(int i = 0) const {return Struct(i)->NoElPerAdcI;}
  float       NoElPerAdcO(int i = 0) const {return Struct(i)->NoElPerAdcO;}
  float       NoElPerAdcX(int i = 0) const {return Struct(i)->NoElPerAdcX;}
  float       OmegaTauScaleI(int i = 0) const {return Struct(i)->OmegaTauScaleI;}
  float       OmegaTauScaleO(int i = 0) const {return Struct(i)->OmegaTauScaleO;}
  float*      SecRowCor(int i = 0)      const {return &Struct(i)->SecRowCorIW[0];}
  float*      SecRowCorIW(int i = 0)    const {return &Struct(i)->SecRowCorIW[0];}
  float*      SecRowCorOW(int i = 0)    const {return &Struct(i)->SecRowCorOW[0];}
  float*      SecRowCorIE(int i = 0)    const {return &Struct(i)->SecRowCorIE[0];}
  float*      SecRowCorOE(int i = 0)    const {return &Struct(i)->SecRowCorOE[0];}

  float*      SecRowSig(int i = 0)      const {return &Struct(i)->SecRowSigIW[0];}
  float*      SecRowSigIW(int i = 0)    const {return &Struct(i)->SecRowSigIW[0];}
  float*      SecRowSigOW(int i = 0)    const {return &Struct(i)->SecRowSigOW[0];}
  float*      SecRowSigIE(int i = 0)    const {return &Struct(i)->SecRowSigIE[0];}
  float*      SecRowSigOE(int i = 0)    const {return &Struct(i)->SecRowSigOE[0];}

  float       PolyaInner(int i = 0)     const {return  Struct(i)->PolyaInner;}
  float       PolyaOuter(int i = 0)     const {return  Struct(i)->PolyaOuter;}
  float       T0offset(int i = 0)       const {return  Struct(i)->T0offset;}
  float       T0offsetI(int i = 0)      const {return  Struct(i)->T0offsetI;}
  float       T0offsetO(int i = 0)      const {return  Struct(i)->T0offsetO;}
  float       FirstRowC(int i = 0)      const {return  Struct(i)->FirstRowC;}
};
#endif
