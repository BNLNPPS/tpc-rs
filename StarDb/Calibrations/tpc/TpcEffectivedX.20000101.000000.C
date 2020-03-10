TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.TpcEffectivedX/TpcEffectivedX Allocated rows: 1  Used rows: 1  Row size: 8 bytes
//  Table: TpcEffectivedX_st[0]--> TpcEffectivedX_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_TpcEffectivedX")) return 0;
TpcEffectivedX_st row;
St_TpcEffectivedX *tableSet = new St_TpcEffectivedX("TpcEffectivedX",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.scaleInner	 =          1; // scale factor for inner dX ;
    row.scaleOuter	 =          1; // -"-        outer dX  ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
