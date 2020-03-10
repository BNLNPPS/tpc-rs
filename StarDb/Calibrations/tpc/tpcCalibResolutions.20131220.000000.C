TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcCalibResolutions/tpcCalibResolutions Allocated rows: 1  Used rows: 1  Row size: 264 bytes
//  Table: tpcCalibResolutions_st[0]--> tpcCalibResolutions_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcCalibResolutions")) return 0;
tpcCalibResolutions_st row;
St_tpcCalibResolutions *tableSet = new St_tpcCalibResolutions("tpcCalibResolutions",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.SpaceCharge	 =       0.05; // SpaceCharge correction ;
    row.GridLeak	 =       -999; // GridLeak correction ;
 memcpy(&row.comment,"First\x20entry\x20from\x20tests\x20with\x20Run\x2014\x20data\x20(entered\x202018-04-09)",60);// comments 
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
