TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcSlowControlSim/tpcSlowControlSim Allocated rows: 1  Used rows: 1  Row size: 112 bytes
//  Table: tpcSlowControlSim_st[0]--> tpcSlowControlSim_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcSlowControlSim")) return 0;
tpcSlowControlSim_st row;
St_tpcSlowControlSim *tableSet = new St_tpcSlowControlSim("tpcSlowControlSim",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.driftVelocity	 =          0; // ;
    row.driftVoltage	 =          0; // ;
    row.innerSectorAnodeVoltage	 =          0; // ;
    row.innerSectorGatingGridV	 =          0; // ;
    row.outerSectorAnodeVoltage	 =          0; // ;
    row.outerSectorGatingGridV	 =          0; // ;
    row.innerSectorGasGain	 =          0; // ;
    row.innerSectorGasGainVzero	 =          0; // ;
    row.innerSectorGasGainb	 =          0; // ;
    row.outerSectorGasGain	 =          0; // ;
    row.outerSectorGasGainVzero	 =          0; // ;
    row.outerSectorGasGainb	 =          0; // ;
    row.hallPressure	 =          0; // ;
    row.hallTemperature	 =          0; // ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
