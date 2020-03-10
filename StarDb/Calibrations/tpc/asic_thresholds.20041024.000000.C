TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.asic_thresholds/asic_thresholds Allocated rows: 1  Used rows: 1  Row size: 16 bytes
//  Table: asic_thresholds_st[0]--> asic_thresholds_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_asic_thresholds")) return 0;
asic_thresholds_st row;
St_asic_thresholds *tableSet = new St_asic_thresholds("asic_thresholds",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.thresh_lo	 =          1; // ;
    row.thresh_hi	 =          4; // ;
    row.n_seq_lo	 =          2; // ;
    row.n_seq_hi	 =          0; // ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
