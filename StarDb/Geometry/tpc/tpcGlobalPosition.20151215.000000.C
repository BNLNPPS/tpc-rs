TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Geometry/tpc/.tpcGlobalPosition/tpcGlobalPosition Allocated rows: 1  Used rows: 1  Row size: 60 bytes
//  Table: tpcGlobalPosition_st[0]--> tpcGlobalPosition_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcGlobalPosition")) return 0;
tpcGlobalPosition_st row;
St_tpcGlobalPosition *tableSet = new St_tpcGlobalPosition("tpcGlobalPosition",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.LocalxShift	 =          0; // cm : x position of TPC center in magnet frame  ;
    row.LocalyShift	 =          0; // cm : y position of TPC center in magnet frame  ;
    row.LocalzShift	 =          0; // cm : z position of TPC center in magnet frame  ;
    row.PhiXY	 =          0; // radians: rotation angle around z axis  (not used) ;
    row.PhiXZ	 =          0; // radians: rotation angle around y axis  XTWIST ;
    row.PhiYZ	 =          0; // radians: rotation angle around x axis  YTWIST ;
    row.XX	 =          1; // XX element of rotation matrix  (not used) ;
    row.YY	 =          1; // YY element of rotation matrix  (not used) ;
    row.ZZ	 =          1; // ZZ element of rotation matrix  (not used) ;
    row.PhiXY_geom	 =          0; // radians: geometrical rotation angle around z axis psi,  -gamma  (not used) ;
    row.PhiXZ_geom	 =          0; // radians: geometrical rotation angle around y axis theta,-beta  ;
    row.PhiYZ_geom	 =          0; // radians: geometrical rotation angle around x axis psi,  -alpha ;
    row.XX_geom	 =          1; // XX element of geometrical rotation matrix  (not used) ;
    row.YY_geom	 =          1; // YY element of geometrical rotation matrix  (not used) ;
    row.ZZ_geom	 =          1; // ZZ element of geometrical rotation matrix  (not used) ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
