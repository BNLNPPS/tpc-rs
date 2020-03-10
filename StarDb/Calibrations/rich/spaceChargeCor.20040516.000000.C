TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/rich/.spaceChargeCor/spaceChargeCor Allocated rows: 1  Used rows: 1  Row size: 64 bytes
//  Table: spaceChargeCor_st[0]--> spaceChargeCor_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_spaceChargeCor")) return 0;
spaceChargeCor_st row;
St_spaceChargeCor *tableSet = new St_spaceChargeCor("spaceChargeCor",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.fullFieldB	 =          0; // Negative Full Field Correction  ;
    row.halfFieldB	 =          0; // Negative Half Field Correction  ;
    row.zeroField	 =          0; // Zero Field " "  ;
    row.halfFieldA	 =          0; // Postive Half " " ;
    row.fullFieldA	 =          0; // Postive Full " " ;
    row.satRate	 =          0; // Saturation Rate Hz  ;
    row.factor	 =          0; // Multiplicative Factor ;
    row.detector	 =          0; // 0=VPDx, 1=BBCx, 2=ZDCx, 3=ZDCe+w, 4=BBCe+w, ... ;
    row.offset	 =          0; // Offset at zero luminosity ;
    row.ewratio	 =          1; // Ratio of charge east/west ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
