TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.tpcPadResponse/tpcPadResponse Allocated rows: 1  Used rows: 1  Row size: 124 bytes
//  Table: tpcPadResponse_st[0]--> tpcPadResponse_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_tpcPadResponse")) return 0;
tpcPadResponse_st row;
St_tpcPadResponse *tableSet = new St_tpcPadResponse("tpcPadResponse",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.innerGasGainFluctuation	 =          0; // unitless  ;
    row.outerGasGainFluctuation	 =          0; // unitless  ;
    row.innerPadResponseSigma	 =          0; // cm  ;
    row.outerPadResponseSigma	 =          0; // cm  ;
    row.innerWirePadCoupling	 =          0; // cm  ;
    row.outerWirePadCoupling	 =          0; // cm  ;
    row.innerRowNormalization	 =          0; // unitless  ;
    row.outerRowNormalization	 =          0; // unitless  ;
    row.BoundaryOfStepFunctions[0]	 =          0; // cm  ;
    row.BoundaryOfStepFunctions[1]	 =          0;
    row.BoundaryOfStepFunctions[2]	 =          0;
    row.BoundaryOfStepFunctions[3]	 =          0;
    row.BoundaryOfStepFunctions[4]	 =          0;
    row.BoundaryOfStepFunctions[5]	 =          0;
    row.innerChargeFractionConstants[0]	 =          0; // unitless  ;
    row.innerChargeFractionConstants[1]	 =          0;
    row.innerChargeFractionConstants[2]	 =          0;
    row.innerChargeFractionConstants[3]	 =          0;
    row.innerChargeFractionConstants[4]	 =          0;
    row.innerChargeFractionConstants[5]	 =          0;
    row.outerChargeFractionConstants[0]	 =          0; // unitless  ;
    row.outerChargeFractionConstants[1]	 =          0;
    row.outerChargeFractionConstants[2]	 =          0;
    row.outerChargeFractionConstants[3]	 =          0;
    row.outerChargeFractionConstants[4]	 =          0;
    row.outerChargeFractionConstants[5]	 =          0;
    row.errorFunctionRange	 =          0; // unitless  ;
    row.errorFunctionEntry	 =          0; // unitless  ;
    row.longitudinalDiffusionConstant	 =          0; // cm/sqrt(cm)  ;
    row.transverseDiffusionConstant	 =          0; // cm/sqrt(cm)  ;
    row.InnerOuterFactor	 =          0; // dimensionless  ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
