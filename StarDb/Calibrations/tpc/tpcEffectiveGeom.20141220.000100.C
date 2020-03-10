TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcEffectiveGeom/tpcEffectiveGeom Allocated rows: 1  Used rows: 1  Row size: 40 bytes
//  Table: tpcEffectiveGeom_st[0]--> tpcEffectiveGeom_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcEffectiveGeom")) return 0;
tpcEffectiveGeom_st row;
St_tpcEffectiveGeom *tableSet = new St_tpcEffectiveGeom("tpcEffectiveGeom",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.drift_length_correction	 =          0; // cm: Diff between actual drift length and  ;
    row.z_inner_offset	 =          0; // cm: Effective distance between  ;
    row.z_outer_offset	 =          0; // cm: Effective distance between  ;
    row.z_inner_offset_West	 =          0; // cm: Effective distance West with respect to East ;
    row.z_outer_offset_West	 =          0; // cm: Effective distance West  -"-                 ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
