/* tpcGas.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
:Description: 
:Synonyms::::
:Source:
:Update:
:Update frequncy:
:Reminder:
:Recall frequency:
:Size of Data:
:Pointer to data:  TPCGas.time:
  type varnam;    //Units : Comments
 mbar : TPC-PT_B
 mbar : TPC-PT_8 @ 7 o'clock
 mbar : TPC-PI_15
 mbar : TPC-PT_7
 degrees K: TPC-T4
 degrees K: TPC-T5
 Gas input values:
 liters/min : TPC-FM_5
 liters/min : TPC-FM_4
 liters/min : TPC-FM_11
 percent    : TPC-CH4_M4
 ppm        : TPC-O2_M1
 Gas exhaust values:
 liters/min : TPC-PT11
 ppm        : TPC-CH4_M3
 ppm        : TPC-H2O_M2
 scf/hr     : TPC-O2_M5
 liters/min : TPC-FI_7
 */
#ifndef TPCGAS_H
#define TPCGAS_H
#define TPCGAS_SPEC \
"struct tpcGas { \
	float barometricPressure; \
	float inputTPCGasPressure; \
	float nitrogenPressure; \
	float gasPressureDiff; \
	float inputGasTemperature; \
	float outputGasTemperature; \
	float flowRateArgon1; \
	float flowRateArgon2; \
	float flowRateMethane; \
	float percentMethaneIn; \
	float ppmOxygenIn; \
	float flowRateExhaust; \
	float percentMethaneOut; \
	float ppmWaterOut; \
	float ppmOxygenOut; \
	float flowRateRecirculation; \
};"
typedef struct tpcGas_st {
	float barometricPressure; 
	float inputTPCGasPressure; 
	float nitrogenPressure; 
	float gasPressureDiff; 
	float inputGasTemperature; 
	float outputGasTemperature; 
	float flowRateArgon1; 
	float flowRateArgon2; 
	float flowRateMethane; 
	float percentMethaneIn; 
	float ppmOxygenIn; 
	float flowRateExhaust; 
	float percentMethaneOut; 
	float ppmWaterOut; 
	float ppmOxygenOut; 
	float flowRateRecirculation; 
} TPCGAS_ST;
#endif /* TPCGAS_H */
