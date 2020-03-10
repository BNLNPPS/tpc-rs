TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcHighVoltages/tpcHighVoltages Allocated rows: 1  Used rows: 1  Row size: 200 bytes
//  Table: tpcHighVoltages_st[0]--> tpcHighVoltages_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcHighVoltages")) return 0;
tpcHighVoltages_st row;
St_tpcHighVoltages *tableSet = new St_tpcHighVoltages("tpcHighVoltages",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.cathode	 =     -27.95; // kVolts  ;
    row.gatedGridRef	 =       -115; // Volts - nominal TPC value but is set by 48 sub-sectors  ;
    row.gridLeakWallTip[0]	 =        111; // Volts - iTPC GridLeak wall tip voltage for 24 sectors  ;
    row.gridLeakWallTip[1]	 =        111;
    row.gridLeakWallTip[2]	 =        111;
    row.gridLeakWallTip[3]	 =        111;
    row.gridLeakWallTip[4]	 =        111;
    row.gridLeakWallTip[5]	 =        111;
    row.gridLeakWallTip[6]	 =        111;
    row.gridLeakWallTip[7]	 =        111;
    row.gridLeakWallTip[8]	 =        111;
    row.gridLeakWallTip[9]	 =        111;
    row.gridLeakWallTip[10]	 =        111;
    row.gridLeakWallTip[11]	 =        111;
    row.gridLeakWallTip[12]	 =        111;
    row.gridLeakWallTip[13]	 =        111;
    row.gridLeakWallTip[14]	 =        111;
    row.gridLeakWallTip[15]	 =        111;
    row.gridLeakWallTip[16]	 =        111;
    row.gridLeakWallTip[17]	 =        111;
    row.gridLeakWallTip[18]	 =        111;
    row.gridLeakWallTip[19]	 =        111;
    row.gridLeakWallTip[20]	 =        111;
    row.gridLeakWallTip[21]	 =        111;
    row.gridLeakWallTip[22]	 =        111;
    row.gridLeakWallTip[23]	 =        111;
    row.gridLeakWallSide[0]	 =        111; // above +100 means no wall  ;
    row.gridLeakWallSide[1]	 =        111;
    row.gridLeakWallSide[2]	 =        111;
    row.gridLeakWallSide[3]	 =        111;
    row.gridLeakWallSide[4]	 =        111;
    row.gridLeakWallSide[5]	 =        111;
    row.gridLeakWallSide[6]	 =        111;
    row.gridLeakWallSide[7]	 =        111;
    row.gridLeakWallSide[8]	 =        111;
    row.gridLeakWallSide[9]	 =        111;
    row.gridLeakWallSide[10]	 =        111;
    row.gridLeakWallSide[11]	 =        111;
    row.gridLeakWallSide[12]	 =        111;
    row.gridLeakWallSide[13]	 =        111;
    row.gridLeakWallSide[14]	 =        111;
    row.gridLeakWallSide[15]	 =        111;
    row.gridLeakWallSide[16]	 =        111;
    row.gridLeakWallSide[17]	 =        111;
    row.gridLeakWallSide[18]	 =        111;
    row.gridLeakWallSide[19]	 =        111;
    row.gridLeakWallSide[20]	 =        111;
    row.gridLeakWallSide[21]	 =        111;
    row.gridLeakWallSide[22]	 =        111;
    row.gridLeakWallSide[23]	 =        111;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
