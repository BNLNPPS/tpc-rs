TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Geometry/tpc/.tpcPadConfig/tpcPadConfig Allocated rows: 1  Used rows: 1  Row size: 24 bytes
//  Table: tpcPadConfig_st[0]--> tpcPadConfig_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcPadConfig")) return 0;
tpcPadConfig_st row;
St_tpcPadConfig *tableSet = new St_tpcPadConfig("tpcPadConfig",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.itpc[0]	 =         0x0; // ;
    row.itpc[1]	 =         0x0;
    row.itpc[2]	 =         0x0;
    row.itpc[3]	 =         0x0;
    row.itpc[4]	 =         0x0;
    row.itpc[5]	 =         0x0;
    row.itpc[6]	 =         0x0;
    row.itpc[7]	 =         0x0;
    row.itpc[8]	 =         0x0;
    row.itpc[9]	 =         0x0;
    row.itpc[10]	 =         0x0;
    row.itpc[11]	 =         0x0;
    row.itpc[12]	 =         0x0;
    row.itpc[13]	 =         0x0;
    row.itpc[14]	 =         0x0;
    row.itpc[15]	 =         0x0;
    row.itpc[16]	 =         0x0;
    row.itpc[17]	 =         0x0;
    row.itpc[18]	 =         0x0;
    row.itpc[19]	 =         0x0;
    row.itpc[20]	 =         0x0;
    row.itpc[21]	 =         0x0;
    row.itpc[22]	 =         0x0;
    row.itpc[23]	 =         0x0;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
