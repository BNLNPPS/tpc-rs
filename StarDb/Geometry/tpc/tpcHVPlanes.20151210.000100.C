TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Geometry/tpc/.tpcHVPlanes/tpcHVPlanes Allocated rows: 1  Used rows: 1  Row size: 36 bytes
//  Table: tpcHVPlanes_st[0]--> tpcHVPlanes_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcHVPlanes")) return 0;
tpcHVPlanes_st row;
St_tpcHVPlanes *tableSet = new St_tpcHVPlanes("tpcHVPlanes",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.CM_shift_z	 =          0; // physical z shift of the CM plane                      ;
    row.CM_tilt_x	 =          0; // x component of the CM plane's normal unit vector      ;
    row.CM_tilt_y	 =          0; // y component of the CM plane's normal unit vector      ;
    row.GGE_shift_z	 =          0; // physical z shift of the GG East plane                 ;
    row.GGE_tilt_x	 =          0; // x component of the GG East plane's normal unit vector ;
    row.GGE_tilt_y	 =          0; // y component of the GG East plane's normal unit vector ;
    row.GGW_shift_z	 =          0; // physical z shift of the GG West plane                 ;
    row.GGW_tilt_x	 =          0; // x component of the GG West plane's normal unit vector ;
    row.GGW_tilt_y	 =          0; // y component of the GG West plane's normal unit vector ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
