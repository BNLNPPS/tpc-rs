TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Conditions/trg/.trgTimeOffset/trgTimeOffset Allocated rows: 1  Used rows: 1  Row size: 12 bytes
//  Table: trgTimeOffset_st[0]--> trgTimeOffset_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_trgTimeOffset")) return 0;
trgTimeOffset_st row;
St_trgTimeOffset *tableSet = new St_trgTimeOffset("trgTimeOffset",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.offset	 =   2.181473; // standard trigger offset in micro-seconds  ;
    row.laserOffset	 =       0.84; // laser trigger offset in micro-seconds  ;
    row.laserOffsetW	 =          0; // laser extra trigger offset for West laser ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
