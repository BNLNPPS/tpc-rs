TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Geometry/tpc/.tpcFieldCage/tpcFieldCage Allocated rows: 1  Used rows: 1  Row size: 12 bytes
//  Table: tpcFieldCage_st[0]--> tpcFieldCage_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcFieldCage")) return 0;
tpcFieldCage_st row;
St_tpcFieldCage *tableSet = new St_tpcFieldCage("tpcFieldCage",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.innerFieldCageShift	 =          0; // cm : z shift of inner field cage w.r.t outer field cage ;
    row.eastClockError	 =          0; // radians :  Phi rotation of East end of TPC in radians ;
    row.westClockError	 =   -0.00043; // radians :  Phi rotation of West end of TPC in radians ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
