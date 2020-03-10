TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcGas/tpcGas Allocated rows: 1  Used rows: 1  Row size: 64 bytes
//  Table: tpcGas_st[0]--> tpcGas_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcGas")) return 0;
tpcGas_st row;
St_tpcGas *tableSet = new St_tpcGas("tpcGas",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.barometricPressure	 =     1018.2; // ;
    row.inputTPCGasPressure	 =       1.99; // ;
    row.nitrogenPressure	 =       1.13; // ;
    row.gasPressureDiff	 =       0.73; // ;
    row.inputGasTemperature	 =     297.66; // ;
    row.outputGasTemperature	 =     297.87; // ;
    row.flowRateArgon1	 =      14.98; // ;
    row.flowRateArgon2	 =       0.44; // ;
    row.flowRateMethane	 =       1.36; // ;
    row.percentMethaneIn	 =      10.17; // ;
    row.ppmOxygenIn	 =      27.39; // ;
    row.flowRateExhaust	 =      11.52; // ;
    row.percentMethaneOut	 =      10.21; // ;
    row.ppmWaterOut	 =       7.86; // ;
    row.ppmOxygenOut	 =      -0.27; // ;
    row.flowRateRecirculation	 =     539.55; // ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
