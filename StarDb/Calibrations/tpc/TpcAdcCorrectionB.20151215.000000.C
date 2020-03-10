TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.TpcAdcCorrectionB/TpcAdcCorrectionB Allocated rows: 2  Used rows: 2  Row size: 120 bytes
//  Table: tpcCorrection_st[0]--> tpcCorrection_st[1]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcCorrection")) return 0;
tpcCorrection_st row;
St_tpcCorrection *tableSet = new St_tpcCorrection("TpcAdcCorrectionB",2);
//
memset(&row,0,tableSet->GetRowSize());
    row.type	 =         11; // ;
    row.idx	 =          1; // ;
    row.nrows	 =          2; // ;
    row.npar	 =          6; // ;
    row.OffSet	 =          0; // ;
    row.min	 =          0; // ;
    row.max	 =          0; // ;
    row.a[0]	 = -0.6673219; // ;
    row.a[1]	 =    45.7099;
    row.a[2]	 = -0.0670727;
    row.a[3]	 =  0.0106947;
    row.a[4]	 =    1.15708;
    row.a[5]	 = -0.0147919;
    row.a[6]	 =          0;
    row.a[7]	 =          0;
    row.a[8]	 =          0;
    row.a[9]	 =          0;
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.type	 =         11; // ;
    row.idx	 =          2; // ;
    row.nrows	 =          2; // ;
    row.npar	 =          6; // ;
    row.OffSet	 =          0; // ;
    row.min	 =          0; // ;
    row.max	 =          0; // ;
    row.a[0]	 = -0.6696278; // ;
    row.a[1]	 =    40.9351;
    row.a[2]	 =  -0.564448;
    row.a[3]	 =   0.132589;
    row.a[4]	 =      1.017;
    row.a[5]	 = 0.00026823;
    row.a[6]	 =          0;
    row.a[7]	 =          0;
    row.a[8]	 =          0;
    row.a[9]	 =          0;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
