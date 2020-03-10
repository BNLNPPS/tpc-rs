TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/rich/.trigDetSums/trigDetSums Allocated rows: 1  Used rows: 1  Row size: 144 bytes
//  Table: trigDetSums_st[0]--> trigDetSums_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_trigDetSums")) return 0;
trigDetSums_st row;
St_trigDetSums *tableSet = new St_trigDetSums("trigDetSums",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.runNumber	 =   16176004; // run number  ;
    row.timeOffset	 = 1435233085; // run begin time  ;
    row.ctbWest	 =          0; // ctb West  ;
    row.ctbEast	 =          0; // ctb East  ;
    row.ctbTOFp	 =          0; // ctbOr + TOFp rate  ;
    row.tofp	 =          0; // TOFp rate  ;
    row.zdcWest	 =          0; // zdc west rate  ;
    row.zdcEast	 =          0; // zdc east rate  ;
    row.zdcX	 =          0; // zdc and rate  ;
    row.mult	 =          0; // mult rate  ;
    row.L0	 =          0; // L0 Rate  ;
    row.bbcX	 =          0; // BBC and Rate  ;
    row.bbcXctbTOFp	 =          0; // BBCAnd + ctbTOFp rate  ;
    row.bbcWest	 =          0; // --BBC West--  ;
    row.bbcEast	 =          0; // --BBC East--  ;
    row.bbcYellowBkg	 =          0; // --(BBC Eastdelayed) and (BBC West)--  ;
    row.bbcBlueBkg	 =          0; // --(BBC Westdelayed) and (BBC East)--  ;
    row.pvpdWest	 =    19127.6; // --PVPD East--  ;
    row.pvpdEast	 =          0; // --PVPD West--  ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
