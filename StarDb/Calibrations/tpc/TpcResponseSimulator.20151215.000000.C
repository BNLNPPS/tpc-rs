TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Calibrations/tpc/.TpcResponseSimulator/TpcResponseSimulator Allocated rows: 1  Used rows: 1  Row size: 212 bytes
//  Table: TpcResponseSimulator_st[0]--> TpcResponseSimulator_st[0]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_TpcResponseSimulator")) return 0;
TpcResponseSimulator_st row;
St_TpcResponseSimulator *tableSet = new St_TpcResponseSimulator("TpcResponseSimulator",1);
//
memset(&row,0,tableSet->GetRowSize());
    row.I0	 =       13.1; // = 13.1 eV, CH4 ;
    row.Cluster	 =        3.2; // = 3.2, average no. of electrons per primary  ;
    row.W	 =       26.2; // = 26.2 eV ;
    row.OmegaTau	 =       3.02; // = 3.02, fit of data ;
    row.K3IP	 =       0.68; // = 0.68,(pads) for a/s = 2.5e-3 and h/s = 0.5 ;
    row.K3IR	 =       0.89; // = 0.89,(row)  for a/s = 2.5e-3 and h/s = 0.5 ;
    row.K3OP	 =       0.55; // = 0.55,(pads) for a/s = 2.5e-3 and h/s = 1.0 ;
    row.K3OR	 =       0.61; // = 0.61,(row)  for a/s = 2.5e-3 and h/s = 1.0 ;
    row.FanoFactor	 =        0.3; // = 0.3 ;
    row.AveragePedestal	 =         50; // = 50.0 ;
    row.AveragePedestalRMS	 =        1.4; // = 1.4, Old Tpc electronics ;
    row.AveragePedestalRMSX	 =        0.7; // = 0.7, New Tpx electronics ;
    row.tauIntegration	 =  1.865e-07; // = 2.5*74.6e-9  secs ;
    row.tauF	 =   3.94e-07; // = 394.0e-9 secs Tpc ;
    row.tauP	 =   7.75e-07; // = 775.0e-9 secs Tpc ;
    row.tauXI	 =      6e-08; // =  60.0e-9 secs Tpx Inner integration time ;
    row.tauXO	 =   7.46e-08; // =  74.6e-9  secs Tpx Outer integration time ;
    row.tauCI	 =          0; // =   0  ;
    row.tauCO	 =          0; // =   0  ;
    row.SigmaJitterTI	 =          0; // = 0.2  for Tpx inner ;
    row.SigmaJitterTO	 =          0; // = 0.2  for Tpx outer ;
    row.SigmaJitterXI	 =          0; // = 0.0  for Tpx inner ;
    row.SigmaJitterXO	 =          0; // = 0.0  for Tpx outer ;
    row.longitudinalDiffusion	 =    0.03624; // cm/sqrt(cm)  ;
    row.transverseDiffusion	 = 0.07056029; // cm/sqrt(cm)  ;
    row.NoElPerAdc	 =        335; // = 335, No. of electrons per 1 ADC count, keep for back compartibility ;
    row.NoElPerAdcI	 =          0; // = 335, No. of electrons per 1 ADC count for inner TPX ;
    row.NoElPerAdcO	 =          0; // = 335, No. of electrons per 1 ADC count for outer TPX ;
    row.NoElPerAdcX	 =          0; // = 335, No. of electrons per 1 ADC count for iTPC      ;
    row.OmegaTauScaleI	 =   3.249675; // ;
    row.OmegaTauScaleO	 =     2.1618; // ;
    row.SecRowCorIW[0]	 =  0.6236582; // parameterization of Inner West correction vs row ;
    row.SecRowCorIW[1]	 = 0.004190659;
    row.SecRowCorOW[0]	 =   1.023561; // parameterization of Outer West correction vs row ;
    row.SecRowCorOW[1]	 = 0.000137881;
    row.SecRowCorIE[0]	 =  0.6236582; // parameterization of Inner East correction vs row ;
    row.SecRowCorIE[1]	 = 0.004190659;
    row.SecRowCorOE[0]	 =   1.023561; // parameterization of Outer East correction vs row ;
    row.SecRowCorOE[1]	 = 0.000137881;
    row.SecRowSigIW[0]	 =  0.0913675; // parameterization of Inner West gain sigma vs row ;
    row.SecRowSigIW[1]	 =          0;
    row.SecRowSigOW[0]	 =  0.0629849; // parameterization of Outer West gain sigma vs row ;
    row.SecRowSigOW[1]	 =          0;
    row.SecRowSigIE[0]	 =  0.0913675; // parameterization of Inner East gain sigma vs row ;
    row.SecRowSigIE[1]	 =          0;
    row.SecRowSigOE[0]	 =  0.0629849; // parameterization of Outer East gain sigma vs row ;
    row.SecRowSigOE[1]	 =          0;
    row.PolyaInner	 =       1.38; // = 1.38, Polya parameter for inner sectors ;
    row.PolyaOuter	 =       1.38; // = 1.38, Polya parameter for outer sectors ;
    row.T0offset	 =   0.356337; // = 0.0   extra off set for Altro chip ;
    row.T0offsetI	 =  0.1258728; // = 0.0   extra off set for inner sector ;
    row.T0offsetO	 = -0.00361778; // = 0.0   extra off set for outer sector ;
    row.FirstRowC	 =          0; // = 0.0   extra correction for the first pad row ;
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
