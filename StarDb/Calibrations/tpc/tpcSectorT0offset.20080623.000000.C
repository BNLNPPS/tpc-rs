TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcSectorT0offset/tpcSectorT0offset Allocated rows: 1  Used rows: 1  Row size: 192 bytes
//  Table: tpcSectorT0offset_st[0]--> tpcSectorT0offset_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcSectorT0offset")) return 0;
tpcSectorT0offset_st row;
St_tpcSectorT0offset *tableSet = new St_tpcSectorT0offset("tpcSectorT0offset",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.t0[0]	 =   -22.2572; // Sector t0 offset per sector: [0-23] Tpx, [24-47] iTpc  (time bins);
    row.t0[1]	 =   -22.2572;
    row.t0[2]	 =   -22.2572;
    row.t0[3]	 =   -22.2572;
    row.t0[4]	 =   -22.2572;
    row.t0[5]	 =   -22.2572;
    row.t0[6]	 =   -22.2572;
    row.t0[7]	 =   -22.2572;
    row.t0[8]	 =   -22.2572;
    row.t0[9]	 =   -22.2572;
    row.t0[10]	 =   -22.2572;
    row.t0[11]	 =   -22.2572;
    row.t0[12]	 =   -22.2572;
    row.t0[13]	 =   -22.2572;
    row.t0[14]	 =   -22.2572;
    row.t0[15]	 =   -22.2572;
    row.t0[16]	 =   -22.2572;
    row.t0[17]	 =   -22.2572;
    row.t0[18]	 =   -22.2572;
    row.t0[19]	 =   -22.2572;
    row.t0[20]	 =   -22.2572;
    row.t0[21]	 =   -22.2572;
    row.t0[22]	 =   -22.2572;
    row.t0[23]	 =   -22.2572;
    row.t0[24]	 =   -22.2572;
    row.t0[25]	 =   -22.2572;
    row.t0[26]	 =   -22.2572;
    row.t0[27]	 =   -22.2572;
    row.t0[28]	 =   -22.2572;
    row.t0[29]	 =   -22.2572;
    row.t0[30]	 =   -22.2572;
    row.t0[31]	 =   -22.2572;
    row.t0[32]	 =   -22.2572;
    row.t0[33]	 =   -22.2572;
    row.t0[34]	 =   -22.2572;
    row.t0[35]	 =   -22.2572;
    row.t0[36]	 =   -22.2572;
    row.t0[37]	 =   -22.2572;
    row.t0[38]	 =   -22.2572;
    row.t0[39]	 =   -22.2572;
    row.t0[40]	 =   -22.2572;
    row.t0[41]	 =   -22.2572;
    row.t0[42]	 =   -22.2572;
    row.t0[43]	 =   -22.2572;
    row.t0[44]	 =   -22.2572;
    row.t0[45]	 =   -22.2572;
    row.t0[46]	 =   -22.2572;
    row.t0[47]	 =   -22.2572;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
