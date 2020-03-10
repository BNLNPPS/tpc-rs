TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcDriftVelocity/tpcDriftVelocity Allocated rows: 1  Used rows: 1  Row size: 16 bytes
//  Table: tpcDriftVelocity_st[0]--> tpcDriftVelocity_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcDriftVelocity")) return 0;
tpcDriftVelocity_st row;
St_tpcDriftVelocity *tableSet = new St_tpcDriftVelocity("tpcDriftVelocity",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.laserDriftVelocityEast	 =    5.56932; // cm/us : from laser beam analysis  ;
    row.laserDriftVelocityWest	 =     5.5679; // cm/us : from laser beam analysis  ;
    row.cathodeDriftVelocityEast	 =          0; // cm/us : from cathode emission  ;
    row.cathodeDriftVelocityWest	 =          0; // cm/us : from cathode emission  ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
