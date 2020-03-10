TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Geometry/tpc/.tpcWirePlanes/tpcWirePlanes Allocated rows: 1  Used rows: 1  Row size: 184 bytes
//  Table: tpcWirePlanes_st[0]--> tpcWirePlanes_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcWirePlanes")) return 0;
tpcWirePlanes_st row;
St_tpcWirePlanes *tableSet = new St_tpcWirePlanes("tpcWirePlanes",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.anodeWireRadius	 =      0.001; // ;
    row.frischGridWireRadius	 =     0.0037; // ;
    row.gatingGridWireRadius	 =     0.0037; // ;
    row.anodeWirePitch	 =        0.4; // ;
    row.frischGridWirePitch	 =        0.1; // ;
    row.gatingGridWirePitch	 =        0.1; // ;
    row.innerSectorAnodeWirePadSep	 =        0.2; // AnodeWire-to-PadPlane distance  ;
    row.innerSectorFrischGridPadSep	 =        0.4; // FrischGrid-to-PadPlane distance  ;
    row.innerSectorGatingGridPadSep	 =          1; // GatingGrid-to-PadPlane distance  ;
    row.outerSectorAnodeWirePadSep	 =        0.4; // AnodeWire-to-PadPlane distance  ;
    row.outerSectorFrischGridPadSep	 =        0.8; // FrischGrid-to-PadPlane distance  ;
    row.outerSectorGatingGridPadSep	 =        1.4; // GatingGrid-to-PadPlane distance  ;
    row.numInnerSectorAnodeWires	 =        170; // ;
    row.numInnerSectorFrischGridWires	 =        681; // ;
    row.numInnerSectorGatingGridWires	 =        681; // ;
    row.firstInnerSectorAnodeWire	 =       53.2; // ;
    row.firstInnerSectorFrischGridWire	 =         53; // ;
    row.firstInnerSectorGatingGridWire	 =         53; // ;
    row.lastInnerSectorAnodeWire	 =      120.8; // ;
    row.numOuterSectorAnodeWires	 =        172; // ;
    row.numOuterSectorFrischGridWires	 =        689; // ;
    row.numOuterSectorGatingGridWires	 =        689; // ;
    row.firstOuterSectorAnodeWire	 =    122.795; // ;
    row.firstOuterSectorFrischGridWire	 =    122.595; // ;
    row.firstOuterSectorGatingGridWire	 =    122.595; // ;
    row.lastOuterSectorAnodeWire	 =    191.195; // ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
