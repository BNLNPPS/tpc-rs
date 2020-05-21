#ifndef tpcrs_structs_h
#define tpcrs_structs_h

struct asic_thresholds {
  int thresh_lo;
  int thresh_hi;
  int n_seq_lo;
  int n_seq_hi;
};

struct g2t_tpc_hit {
  int id; /* primary key */
  int next_tr_hit_p; /* Id of next hit on same track */
  int track_p; /* Id of parent track */
  int volume_id; /* STAR volume identification */
  float de; /* energy deposition at hit */
  float ds; /* path length within padrow */
  float p[3]; /* local momentum */
  float tof; /* time of flight */
  float x[3]; /* coordinate (Cartesian) */
  float lgam; /* ALOG10(GEKin/AMass) */
  float length; /* track length up to this hit */
  float adc; /* signal in ADC  after digitization */
  float pad; /* hit pad position used in digitization */
  float timebucket; /* hit time position -"- */
  int np; /* no. of primary electrons */
};

struct g2t_track {
  int eg_label; /* generator track label (0 if GEANT track) */
  int eg_pid; /* event generator particle id */
  int ge_pid; /* GEANT particle id */
  int hit_ctb_p; /* Id of first ctb hit on track linked list */
  int hit_eem_p; /* Id of first eem hit on track linked list */
  int hit_emc_p; /* Id of first emc hit on track linked list */
  int hit_esm_p; /* Id of first esm hit on track linked list */
  int hit_ftp_p; /* Id of first ftp hit on track linked list */
  int hit_gem_p; /* Id of first gem hit on track linked list */
  int hit_hpd_p; /* Id of first hpd hit on track linked list */
  int hit_ist_p; /* Id of first ist hit on track linked list */
  int hit_igt_p; /* Id of first igt hit on track linked list */
  int hit_fst_p; /* Id of first fst hit on track linked list */
  int hit_fgt_p; /* Id of first fgt hit on track linked list */
  int hit_fpd_p; /* Id of first fpd hit on track linked list */
  int hit_fsc_p; /* Id of first fsc hit on track linked list */
  int hit_mtd_p; /* Id of first mtd hit on track linked list */
  int hit_mwc_p; /* Id of first mwc hit on track linked list */
  int hit_pgc_p; /* Id of first pgc hit on track linked list */
  int hit_pmd_p; /* Id of first psc hit on track linked list */
  int hit_smd_p; /* pointer to first SHM linked list hit */
  int hit_ssd_p; /* Id of first ssd hit on track linked list */
  int hit_svt_p; /* Id of first svt hit on track linked list */
  int hit_pix_p; /* Id of first pix hit on track linked list */
  int hit_tof_p; /* Id of first tof hit on track linked list */
  int hit_tpc_p; /* Id of first tpc hit on track linked list */
  int hit_vpd_p; /* Id of first vpd hit on track linked list */
  int hit_etr_p; /* Id of first etr hit on track linked list */
  int hit_hca_p; /* Id of first hca hit on track linked list */
  int hit_fts_p; /* Id of first fts hit on track linked list */
  int hit_eto_p; /* Id of first etof hit on track linked list */
  int id; /* primary key */
  int is_shower; /* 1 if shower track, 0 if not */
  int itrmd_vertex_p; /* First intermediate vertex */
  int n_ctb_hit; /* Nhits in ctb */
  int n_eem_hit; /* Nhits in eem (endcap em cal) */
  int n_emc_hit; /* Nhits in emc */
  int n_esm_hit; /* Nhits in esm (endcap shower max) */
  int n_ftp_hit; /* Nhits in forward tpc */
  int n_gem_hit; /* Nhits in gem barrel */
  int n_hpd_hit; /* Nhits in hpd */
  int n_ist_hit; /* Nhits in ist */
  int n_igt_hit; /* Nhits in igt */
  int n_fst_hit; /* Nhits in fst [new Forward Silicon Tracker] */
  int n_fgt_hit; /* Nhits in fgt */
  int n_fpd_hit; /* Nhits in fpd */
  int n_fsc_hit; /* Nhits in fsc */
  int n_mtd_hit; /* Nhits in mtd */
  int n_mwc_hit; /* Nhits in mwc */
  int n_pgc_hit; /* Nhits in pgc  ???  */
  int n_pmd_hit; /* Nhits in pmd (PMD) */
  int n_smd_hit; /* number of hits in shower max */
  int n_ssd_hit; /* Nhits in ssd */
  int n_svt_hit; /* Nhits in svt */
  int n_pix_hit; /* Nhits in pix */
  int n_tof_hit; /* Nhits in tof */
  int n_tpc_hit; /* Nhits in tpc + (no. of real hits after TpcRS) << 8 */
  int n_vpd_hit; /* Nhits in vpd */
  int n_etr_hit; /* Nhits on etr */
  int n_hca_hit; /* Nhits on hca [hadronic calorimeter]*/
  int n_fts_hit; /* Nhits on fts */
  int n_eto_hit; /* Nhits on etof */
  int n_stg_hit; /* Nhits on sTGC [new Forward small thin gap chambers]*/
  int n_wca_hit; /* Nhits on wca [West EM calorimeter]*/
  int next_parent_p; /* Id of next parent track */
  int next_vtx_trk_p; /* Next daughter track of start vertex */
  int start_vertex_p; /* Id of start vertex of track */
  int stop_vertex_p; /* Id of stop vertex of this track */
  float charge; /* Charge */
  float e; /* Energy */
  float eta; /* Pseudorapidity */
  float p[3]; /* Momentum */
  float pt; /* Transverse momentum */
  float ptot; /* Total momentum */
  float rapidity; /* Rapidity */
};

struct g2t_vertex {
  char ge_volume[4]; /* GEANT volume name */
  int daughter_p; /* Id of first daughter in linked list */
  int eg_label; /* generator label (or 0 if GEANT vertex) */
  int eg_proc; /* generator production mechanism (if any) */
  int event_p; /* pointer to event */
  int ge_medium; /* GEANT Medium */
  int ge_proc; /* >0 GEANT, =0 event generator */
  int id; /* primary key */
  int is_itrmd; /* flags intermediate vertex */
  int n_daughter; /* Number of daughter tracks */
  int n_parent; /* number of parent tracks */
  int next_itrmd_p; /* Id of next intermedate vertex */
  int next_prim_v_p; /* Id of next primary vertex */
  int parent_p; /* Id of first parent track */
  float eg_tof; /* generator vertex production time */
  float eg_x[3]; /* generator vertex coordinate (Cartesian) */
  float ge_tof; /* GEANT vertex production time */
  float ge_x[3]; /* GEANT vertex coordinate (Cartesian) */
};

struct itpcPadGainT0 {
  int run; /* pulser run number used */
  float Gain[24][40][120]; /* Gains per pad*/
  float T0[24][40][120]; /* T0 per pad*/
};

/** Survey data */
struct iTPCSurvey {
	int Id; 
	float Angle; 
	float dx; 
	float dy; 
	float ScaleX; 
	float ScaleY; 
	char comment[32]; 
};

struct MagFactor {
	float ScaleFactor; 
};

/** TMultiDimFit
 row index
 total no. of real rows in the table; For Db interface (where nrows = 50)
 type = 0 kMonomials, type = 1 kChebyshev, type = 2 kLegendre
 == 2 for now.
  p_ij = Power[i * NVariables + j];
 */
struct MDFCorrection {
	unsigned char idx; 
	unsigned char nrows; 
	unsigned char PolyType; 
	unsigned char NVariables; 
	unsigned char NCoefficients; 
	unsigned char Power[100]; 
	double DMean; 
	double XMin[2]; 
	double XMax[2]; 
	double Coefficients[50]; 
	double CoefficientsRMS[50]; 
};

struct richvoltages {
  unsigned int runNumber; /*   */
  unsigned int startStatusTime; /*   */
  unsigned int endStatusTime; /*   */
  unsigned int status; /*   */
};

/** Table for Space Charge Corrections */
struct spaceChargeCor {
	double fullFieldB; /* Negative Full Field Correction  */
	double halfFieldB; /* Negative Half Field Correction  */
	double zeroField; /*  Zero Field " "  */
	double halfFieldA; /*  Postive Half " " */
	double fullFieldA; /*  Postive Full " " */
	double satRate; /* Saturation Rate Hz  */
	float factor; /*  Multiplicative Factor */
	float detector; /* 0=VPDx, 1=BBCx, 2=ZDCx, 3=ZDCe+w, 4=BBCe+w, ... */
	float offset; /* Offset at zero luminosity */
	float ewratio; /* Ratio of charge east/west */
};

struct starClockOnl {
  unsigned int runNumber; /*   run number  */
  unsigned int time; /*   unix time of entry  */
  double frequency; /*   frequency in Hz  */
};

struct starMagOnl {
	unsigned int runNumber; /*   run number  */
	unsigned int time; /*   unix time of entry  */
	double current; /*   magnet current (- means B polarity)  */
};

/** Survey data
   Translation from local to Master
   m = R*l + t
           ( r00 r01 r02 ) (xl)   ( t0 )    ( r00 r01 r02 ) (xl)   ( t0 )    (     1 -gamma   beta ) (xl)   ( t0 )
       R = ( r10 r11 r12 ) (yl) + ( t1 ) =  ( r10 r11 r12 ) (zl) + ( t1 ) ~  ( gamma      1 -alpha ) (zl) + ( t1 ) 
           ( r20 r21 r22 ) (zl)   ( t2 )    ( r20 r21 r22 ) (yl)   ( t2 )    ( -beta  alpha      1 ) (yl)   ( t2 )
 SVT
 Id = 0                                for SvtOnGlobal
 Id = [0, 1]                           for ShellOnGlobal, 
 0 is the x (South) Shell, 1 is the -x (North) Shell"
 Id = 1000*barrel + ladder             for LadderOnSurvey
 Id = 1000*barrel + ladder             for LadderOnShell
 Id = 1000*barrel + 100*wafer + ladder for WaferOnLadder
 SSD
 Id = 0                                for SsdOnGlobal
 Id = sector [1-4]                     SsdSectorsOnGlobal
 Id = 100*sector + ladder              SsdLaddersOnSectors
 Id = 7000 + 100*wafer + ladder        SsdWafersOnLadders
 */
struct Survey {
	int Id; 
	double r00; 
	double r01; /* -gamma */
	double r02; /*  beta  */
	double r10; /*  gamma */
	double r11; 
	double r12; /* -alpha */
	double r20; /* -beta  */
	double r21; /*  alpha */
	double r22; 
	double t0; 
	double t1; 
	double t2; 
	double sigmaRotX; 
	double sigmaRotY; 
	double sigmaRotZ; 
	double sigmaTrX; 
	double sigmaTrY; 
	double sigmaTrZ; 
	char comment[32]; 
};

/** Parameters of Altro cheap configuration
        N  no. of Altro parameters,
        N  0  > filter is switched off
        N  1 > old TPC electronics
        Altro::ConfigZerosuppression(int Threshold, int MinSamplesaboveThreshold, int Presamples, int Postsamples)
        Altro::ConfigTailCancellationFilter(int K1, int K2, int K3, int L1, int L2, int L3)
  Threshold
  MinSamplesaboveThreshold
  K1 coefficient of the TCF
  K2 coefficient of the TCF
  K3 coefficient of the TCF
  L1 coefficient of the TCF
  L2 coefficient of the TCF
  L3 coefficient of the TCF
 */
struct tpcAltroParams {
  int N; /*  = no. of Altro parameters */
  int Altro_thr;
  int Altro_seq;
  int Altro_K1;
  int Altro_K2;
  int Altro_K3;
  int Altro_L1;
  int Altro_L2;
  int Altro_L3;
};

/** average voltages over stable periods of time, */
struct tpcAnodeHVavg {
  unsigned short sector; /*  sector 1-24 */
  unsigned short socket; /*  MWC socket/card (ISOR=17,OSIR=18,OSOR=19)  */
  float voltage; /*  average voltage  */
  float rms; /*  rms for averaged voltage */
  int numentries; /*  number of entries used for average */
  int numoutliers; /*  number of encountered outliers */
};

struct tpcAnodeHV {
  unsigned short sector; /*  sector 1-24 */
  unsigned short socket; /*  MWC socket/card (ISOR=17,OSIR=18,OSOR=19)  */
  float voltage; /*   HV setting  */
};

struct TpcAvgCurrent {
  int run; /* run no. used for averaging  */
  int start_time; /* begin unix time of averaging interval */
  int stop_time; /* end   unix time of averaging interval */
  float AvCurrent[192]; /* average current per sector(24) and channel(8) [muA]*/
  float AcCharge[192]; /* accumulated charge per sector(24) and channel(8) [C]*/
};

struct TpcAvgPowerSupply {
  int run; /* run no. used for averaging  */
  int start_time; /* begin unix time of averaging interval */
  int stop_time; /* end   unix time of averaging interval */
  float Current[192]; /* average current per sector(24) and channel(8) [muA]*/
  float Charge[192]; /* accumulated charge per sector(24) and channel(8) [C]*/
  float Voltage[192]; /* average Voltage per sector(24) and channel(8) [V]*/
};

/** Table for Resolutions of TPC Calibrations */
struct tpcCalibResolutions {
	float SpaceCharge; /* SpaceCharge correction */
	float GridLeak; /* GridLeak correction */
	char comment[255]; /* comments */
};

/** Table of events depositing significant charge into the TPC;
                    bunchCrossing is stored as 32bit pairs because no 64bit integer
                    works with the database, i.e. the first bunch crissing is actually
                    (eventBunchCrossingsHigh[0] << 32)  eventBunchCrossingsLow[0]
                    with appropriate 32to64 conversion before bitshifting
 */
struct tpcChargeEvent {
	int nChargeEvents; /* number of charge events in this record */
	unsigned int eventBunchCrossingsLow[4096]; /* number of bunches into the run when charge event occurred (32 low bits) */
	unsigned int eventBunchCrossingsHigh[4096]; /* number of bunches into the run when charge event occurred (32 high bits) */
	float eventCharges[4096]; /* metric of magnitude of charge deposited in TPC */
	int badBunch; /* collider bunch where most of the charge events occurred */
};

/** Drift Distance depended correction
 type = 0 polymonical fit,                                        use only [min,max]
 type = 1 TChebyshev poly in range [min,max] => [-1,1]
 type = 2 shifted  TChebyshev poly in range [min,max] => [ 0,1]
 type = 3 X => Log(1 - |x|)                                       use only [min,max]
 type = 4 X => Log(1 - |x|)*sign(x)                                -"-
 type = 5 X => Log(x), for x <= 1 => 0
 type = 10 (== log(1. + OffSet/x) + poly(x,npar))                  -"-
 type = 11 (== log(1. + OffSet/x) + poly(x,npar) for log(ADC) and |Z|
 type = 200 cut on range [min,max]
 type = 300 don't correct out of range [min,max]
 type = 1000      ; gaus(0)+pol0(3);
 type = 1000 + 100; gaus(0)+pol1(3)
 type = 1000 + 200; gaus(0)+pol2(3)
 type = 1000 + 300; gaus(0)+pol3(3)
 type = 2000      ; expo(0)+pol0(2);
 type = 2000 + 100; expo(0)+pol1(2)
 type = 2000 + 200; expo(0)+pol2(2)
 type = 2000 + 300; expo(0)+pol3(2)
 row index
COMMENTS TRUNCATED */
struct tpcCorrection {
  int type;
  int idx;
  int nrows;
  int npar;
  double OffSet;
  double min;
  double max;
  double a[10];
};

struct tpcDimensions {
  int numberOfSectors; /*   */
  double tpcInnerRadius; /*   */
  double tpcOuterRadius; /*   */
  double tpcTotalLength; /*   */
  double wheelInnerRadius; /*   */
  double wheelOuterRadius; /*   */
  double wheelThickness; /*   */
  double senseGasOuterRadius; /*   */
  double tpeaThickness; /*   */
  double cathodeInnerRadius; /*   */
  double cathodeOuterRadius; /*   */
  double cathodeThickness; /*   */
  double outerCuThickness; /*   */
  double outerKaptonThickness; /*   */
  double outerNomexThickness; /*   */
  double outerGlueThickness; /*   */
  double outerInsGasThickness; /*   */
  double outerAlThickness; /*   */
  double outerAlHoneycombThickness; /*   */
  double innerGlueThickness; /*   */
  double innerNomexThickness; /*   */
  double innerKaptonThickness; /*   */
  double innerAlThickness; /*   */
  double innerGapWidI; /*  inner width of air in support wheel  */
  double innerGapWidO; /*  outer width of air in support wheel  */
  double innerGapHeit; /*  height (dr) of air in support wheel  */
  double innerGapRad; /*  air in support wheel - center radius  */
  double innerInWidth; /*   sector width at inner radius  */
  double innerOutWidth; /*   sector width at outer radius  */
  double innerHeight; /*   sector radial height  */
  double innerPPDepth; /*  padplane thickness (Al+pcb)  */
  double innerAlDepth; /*  depth of openings in Al structure  */
  double innerMWCDepth; /*  full thickness MWC sensitive region  */
  double innerBoundary; /*  Al frame boundary width  */
  double innerRCenter; /*  sector center radius (precision hole)  */
  double innerMWCInn; /*  MWC sensitive region inner size  */
  double innerMWCOut; /*  MWC sensitive region outer size  */
  double innerMVCHei; /*  MVC sensitive region radial  */
  int innerAirGaps; /*   number of air gaps in Al  */
  int innerExtraAl; /*   number of extra Al supportpieces  */
  double innerZGaps[5]; /*   opening positions  */
  double innerZGapsSize[5]; /*   opening size  */
  double innerXExtraAl[5]; /*   x positions  */
  double innerZExtraAl[5]; /*   z positions  */
  double innerDXExtraAl[5]; /*   x thickness  */
  double innerDZExtraAl[5]; /*   z thickness  */
  double outerGapWidI; /*  inner width of air in support wheel  */
  double outerGapWidO; /*  outer width of air in support wheel  */
  double outerGapHeit; /*  height (dr) of air in support wheel  */
  double outerGapRad; /*  air in support wheel - center radius  */
  double outerInWidth; /*   sector width at inner radius  */
  double outerOutWidth; /*   sector width at outer radius  */
  double outerHeight; /*   sector radial height  */
  double outerPPDepth; /*  padplane thickness (Al+pcb)  */
  double outerAlDepth; /*  depth of openings in Al structure  */
  double outerMWCDepth; /*  full thickness MWC sensitive region  */
  double outerBoundary; /*  Al frame boundary width  */
  double outerRCenter; /*  sector center radius (precision hole)  */
  double outerMWCInn; /*  MWC sensitive region inner size  */
  double outerMWCOut; /*  MWC sensitive region outer size  */
  double outerMVCHei; /*  MVC sensitive region radial  */
  int outerAirGaps; /*   number of air gaps in Al  */
  int outerExtraAl; /*   number of extra Al supportpieces  */
  double outerZGaps[8]; /*   opening positions  */
  double outerZGapsSize[8]; /*   opening size  */
  double outerXExtraAl[5]; /*   x positions  */
  double outerZExtraAl[5]; /*   z positions  */
  double outerDXExtraAl[5]; /*   x thickness  */
  double outerDZExtraAl[5]; /*   z thickness  */
};

struct tpcDriftVelocity {
	float laserDriftVelocityEast; /*   cm/us : from laser beam analysis  */
	float laserDriftVelocityWest; /*   cm/us : from laser beam analysis  */
	float cathodeDriftVelocityEast; /*   cm/us : from cathode emission  */
	float cathodeDriftVelocityWest; /*   cm/us : from cathode emission  */
};

/** Effective height of pad row */
struct TpcEffectivedX {
	float scaleInner; /* scale factor for inner dX */
	float scaleOuter; /*       -"-        outer dX  */
};

struct tpcEffectiveGeom {
  double drift_length_correction; /*  cm: Diff between actual drift length and  */
  double z_inner_offset; /*  cm: Effective distance between  */
  double z_outer_offset; /*  cm: Effective distance between  */
  double z_inner_offset_West; /*  cm: Effective distance West with respect to East */
  double z_outer_offset_West; /*  cm: Effective distance West  -"-                 */
};

struct tpcElectronics {
  int numberOfTimeBins; /*   */
  double nominalGain; /*   mV/fC  */
  double samplingFrequency; /* MHz, not used,  overwritten by starClockOnl*/
  double tZero; /*   us (microseconds)  */
  double adcCharge; /*   fC/adc count  */
  double adcConversion; /*   mV/adc count  */
  double averagePedestal; /*   adc counts  */
  double shapingTime; /*   ns  */
  double tau; /*   ns  */
};

struct tpcFieldCage {
  float innerFieldCageShift; /* cm : z shift of inner field cage w.r.t outer field cage */
  float eastClockError; /* radians :  Phi rotation of East end of TPC in radians */
  float westClockError; /* radians :  Phi rotation of West end of TPC in radians */
};

/** provide info on shorted rings in the field cages */
struct tpcFieldCageShort {
	float side; /* 0 = east, 1 = west */
	float cage; /* 0 = inner, 1 = outer */
	float ring; /* ring location of the short (e.g. 169.5) */
	float resistor; /* MOhm value of added external resistor to resistor chain */
	float MissingResistance; /* missing resistance */
};

/**
  type varnam;    //Units : Comments
 mbar : TPC-PT_B
 mbar : TPC-PT_8 @ 7 o'clock
 mbar : TPC-PI_15
 mbar : TPC-PT_7
 degrees K: TPC-T4
 degrees K: TPC-T5
 Gas input values:
 liters/min : TPC-FM_5
 liters/min : TPC-FM_4
 liters/min : TPC-FM_11
 percent    : TPC-CH4_M4
 ppm        : TPC-O2_M1
 Gas exhaust values:
 liters/min : TPC-PT11
 ppm        : TPC-CH4_M3
 ppm        : TPC-H2O_M2
 scf/hr     : TPC-O2_M5
 liters/min : TPC-FI_7
 */
struct tpcGas {
  float barometricPressure;
  float inputTPCGasPressure;
  float nitrogenPressure;
  float gasPressureDiff;
  float inputGasTemperature;
  float outputGasTemperature;
  float flowRateArgon1;
  float flowRateArgon2;
  float flowRateMethane;
  float percentMethaneIn;
  float ppmOxygenIn;
  float flowRateExhaust;
  float percentMethaneOut;
  float ppmWaterOut;
  float ppmOxygenOut;
  float flowRateRecirculation;
};

struct tpcGlobalPosition {
  float LocalxShift; /* cm : x position of TPC center in magnet frame  */
  float LocalyShift; /* cm : y position of TPC center in magnet frame  */
  float LocalzShift; /* cm : z position of TPC center in magnet frame  */
  float PhiXY; /* radians: rotation angle around z axis  (not used) */
  float PhiXZ; /* radians: rotation angle around y axis  XTWIST */
  float PhiYZ; /* radians: rotation angle around x axis  YTWIST */
  float XX; /* XX element of rotation matrix  (not used) */
  float YY; /* YY element of rotation matrix  (not used) */
  float ZZ; /* ZZ element of rotation matrix  (not used) */
  float PhiXY_geom; /* radians: geometrical rotation angle around z axis psi,  -gamma  (not used) */
  float PhiXZ_geom; /* radians: geometrical rotation angle around y axis theta,-beta  */
  float PhiYZ_geom; /* radians: geometrical rotation angle around x axis psi,  -alpha */
  float XX_geom; /* XX element of geometrical rotation matrix  (not used) */
  float YY_geom; /* YY element of geometrical rotation matrix  (not used) */
  float ZZ_geom; /* ZZ element of geometrical rotation matrix  (not used) */
};

/** Table for TPC GridLeaks */
struct tpcGridLeak {
	double InnerGLRadius; /* Radius of GL around inner sectors           */
	double MiddlGLRadius; /* Radius of GL between inner/outer sectors    */
	double OuterGLRadius; /* Radius of GL around outer sectors           */
	double InnerGLWidth; /* Width of GL around inner sectors            */
	double MiddlGLWidth; /* Width of GL between inner/outer sectors     */
	double OuterGLWidth; /* Width of GL around outer sectors            */
	double InnerGLStrength; /* Strength of GL around inner sectors         */
	double MiddlGLStrength; /* Strength of GL between inner/outer sectors  */
	double OuterGLStrength; /* Strength of GL around outer sectors         */
};

/** Cathode and gating grid voltage for ExB distortions  Cathode and gating grid voltage for ExB distortions */
struct tpcHighVoltages {
	float cathode; /*   kVolts  */
	float gatedGridRef; /*   Volts - nominal TPC value but is set by 48 sub-sectors  */
	float gridLeakWallTip[24]; /*   Volts - iTPC GridLeak wall tip voltage for 24 sectors  */
	float gridLeakWallSide[24]; /*   above +100 means no wall  */
};

/** orientation of TPC HV planes:
                    central membrane
                    east gated grid
                    west gated grid
 */
struct tpcHVPlanes {
	float CM_shift_z; /* physical z shift of the CM plane                      */
	float CM_tilt_x; /* x component of the CM plane's normal unit vector      */
	float CM_tilt_y; /* y component of the CM plane's normal unit vector      */
	float GGE_shift_z; /* physical z shift of the GG East plane                 */
	float GGE_tilt_x; /* x component of the GG East plane's normal unit vector */
	float GGE_tilt_y; /* y component of the GG East plane's normal unit vector */
	float GGW_shift_z; /* physical z shift of the GG West plane                 */
	float GGW_tilt_x; /* x component of the GG West plane's normal unit vector */
	float GGW_tilt_y; /* y component of the GG West plane's normal unit vector */
};

/** Table for TPC OmegaTau */
struct tpcOmegaTau {
	float tensorV1; /* tensor for OmegaTau           */
	float tensorV2; /* tensor for OmegaTau    */
	unsigned short distortionCorrectionsMode; /* modes for field distortion and non-uniformity corrections  */
};

/**
TPC padrow configurations (tpcitpc) for 24 sectors
 sector status 0 => tpc, 1 => itpc
 */
struct tpcPadConfig {
  unsigned char itpc[24];
};

struct tpcPadGainT0 {
  int run; /* pulser run number used */
  float Gain[24][45][182]; /* Gains per pad*/
  float T0[24][45][182]; /* T9 per pad*/
};

struct tpcPadPlanes {
  int padRows; /*   */
  int innerPadRows; /*   */
  int innerPadRows48; /*   */
  int innerPadRows52; /*   */
  int outerPadRows; /*   */
  int superInnerPadRows; /*   */
  int superOuterPadRows; /*   */
  double innerSectorPadWidth; /*   */
  double innerSectorPadLength; /*   */
  double innerSectorPadPitch; /*   */
  double innerSectorRowPitch1; /*   */
  double innerSectorRowPitch2; /*   */
  double firstPadRow; /*   */
  double firstOuterSectorPadRow; /*   */
  double lastOuterSectorPadRow; /*   */
  double firstRowWidth; /*   */
  double lastRowWidth; /*   */
  double outerSectorPadWidth; /*   */
  double outerSectorPadLength; /*   */
  double outerSectorPadPitch; /*   */
  double outerSectorRowPitch; /*   */
  double outerSectorLength; /*   */
  double ioSectorSeparation; /*   */
  double innerSectorEdge; /*   */
  double outerSectorEdge; /*   */
  double innerSectorPadPlaneZ; /*   */
  double outerSectorPadPlaneZ; /*   */
  int innerPadsPerRow[13]; /*   */
  int outerPadsPerRow[32]; /*   */
  double innerRowRadii[13]; /*   */
  double outerRowRadii[32]; /*   */
};

/** T0s by padrow, applied after clusterfinding */
struct tpcPadrowT0 {
	float T0[100]; /*  T0s per padrow */
};

/**
Map tpc row, pad range [padMin, padMax] to rdo number
 row index
 total no. of real rows in the table; For Db interface (where nrows = 50)
 */
struct tpcRDOMap {
	int idx; 
	int nrows; 
	unsigned char row; 
	unsigned char padMin; 
	unsigned char padMax; 
	unsigned char rdo; 
};

struct tpcRDOMasks {
  unsigned int runNumber; /*       run number  */
  unsigned int sector; /*   sector  */
  unsigned int mask; /*   enable mask  */
};

struct tpcRDOT0offset {
	unsigned char isShifted[24]; /* flag if there is any RDO off set to tsector */
	float t0[24][10]; /* RDO t0 offset per sector: [0-23], rdo [0-5] Tpx, [6-9] iTpc  (time bins) */
};

/**
Tpc Response Simulator parameters  Tpc Response Simulator parameters
  2.145, effective reduction of OmegaTau near
				  Inner sector anode wire
  1.8, effective reduction of OmegaTau near
				  Outer sector anode wire
 */
struct TpcResponseSimulator {
  float I0; /* = 13.1 eV, CH4 */
  float Cluster; /* = 3.2, average no. of electrons per primary  */
  float W; /* = 26.2 eV */
  float OmegaTau; /* = 3.02, fit of data */
  float K3IP; /* = 0.68,(pads) for a/s = 2.5e-3 and h/s = 0.5 */
  float K3IR; /* = 0.89,(row)  for a/s = 2.5e-3 and h/s = 0.5 */
  float K3OP; /* = 0.55,(pads) for a/s = 2.5e-3 and h/s = 1.0 */
  float K3OR; /* = 0.61,(row)  for a/s = 2.5e-3 and h/s = 1.0 */
  float FanoFactor; /* = 0.3 */
  float AveragePedestal; /* = 50.0 */
  float AveragePedestalRMS; /* = 1.4, Old Tpc electronics */
  float AveragePedestalRMSX; /* = 0.7, New Tpx electronics */
  float tauIntegration; /* = 2.5*74.6e-9  secs */
  float tauF; /* = 394.0e-9 secs Tpc */
  float tauP; /* = 775.0e-9 secs Tpc */
  float tauXI; /* =  60.0e-9 secs Tpx Inner integration time */
  float tauXO; /* =  74.6e-9  secs Tpx Outer integration time */
  float tauCI; /* =   0  */
  float tauCO; /* =   0  */
  float SigmaJitterTI; /* = 0.2  for Tpx inner */
  float SigmaJitterTO; /* = 0.2  for Tpx outer */
  float SigmaJitterXI; /* = 0.0  for Tpx inner */
  float SigmaJitterXO; /* = 0.0  for Tpx outer */
  float longitudinalDiffusion; /*   cm/sqrt(cm)  */
  float transverseDiffusion; /*   cm/sqrt(cm)  */
  float NoElPerAdc; /* = 335, No. of electrons per 1 ADC count, keep for back compartibility */
  float NoElPerAdcI; /* = 335, No. of electrons per 1 ADC count for inner TPX */
  float NoElPerAdcO; /* = 335, No. of electrons per 1 ADC count for outer TPX */
  float NoElPerAdcX; /* = 335, No. of electrons per 1 ADC count for iTPC      */
  float OmegaTauScaleI;
  float OmegaTauScaleO;
  float SecRowCorIW[2]; /* parameterization of Inner West correction vs row */
  float SecRowCorOW[2]; /* parameterization of Outer West correction vs row */
  float SecRowCorIE[2]; /* parameterization of Inner East correction vs row */
  float SecRowCorOE[2]; /* parameterization of Outer East correction vs row */
  float SecRowSigIW[2]; /* parameterization of Inner West gain sigma vs row */
  float SecRowSigOW[2]; /* parameterization of Outer West gain sigma vs row */
  float SecRowSigIE[2]; /* parameterization of Inner East gain sigma vs row */
  float SecRowSigOE[2]; /* parameterization of Outer East gain sigma vs row */
  float PolyaInner; /* = 1.38, Polya parameter for inner sectors */
  float PolyaOuter; /* = 1.38, Polya parameter for outer sectors */
  float T0offset; /* = 0.0   extra off set for Altro chip */
  float T0offsetI; /* = 0.0   extra off set for inner sector */
  float T0offsetO; /* = 0.0   extra off set for outer sector */
  float FirstRowC; /* = 0.0   extra correction for the first pad row */
};

/**
Table for SpaceCharge and GridLeak Correction parameters
       SC parameters: 4 array elements for west, plus 4 for east  8,
          each of the 4 elements are additive terms in the SC formula
       GL parameters: 24 array elements for TPC sectors,
          plus radius and width of the charge sheet
       scaler and mode definitions in StDetectorDbMakerSt_tpcSCGLC.h
 */
struct tpcSCGL {
	float SC[8]; /* Scale factor relating luminosity scaler to SpaceCharge */
	float SCoffset[8]; /* Offset to define luminosity for SpaceCharge */
	float SCexponent[8]; /* Luminosity exponential factor for SpaceCharge */
	float SCscaler[8]; /* Luminosity detector scaler */
	float GL[24]; /* Scale factor relating SpaceCharge to GridLeak */
	float GLoffset[24]; /* Offset to define luminosity for GridLeak */
	float GLradius; /* Radius of GridLeak between inner/outer sectors */
	float GLwidth; /* Width of GridLeak between inner/outer sectors */
	int mode; /* Modes to simplify parameter controls */
	char comment[256]; 
};

/** Tpc gain correction for sector  row */
struct TpcSecRowCor {
	float GainScale[100]; /*  Gains for sector & row */
	float GainRms[100]; /*  RMS  - " -  */
};

struct tpcSectorT0offset {
	float t0[48]; /* Sector t0 offset per sector: [0-23] Tpx, [24-47] iTpc  (time bins)*/
};

struct tpcWirePlanes {
  double anodeWireRadius; /*   */
  double frischGridWireRadius; /*   */
  double gatingGridWireRadius; /*   */
  double anodeWirePitch; /*   */
  double frischGridWirePitch; /*   */
  double gatingGridWirePitch; /*   */
  double innerSectorAnodeWirePadSep; /*   AnodeWire-to-PadPlane distance  */
  double innerSectorFrischGridPadSep; /*   FrischGrid-to-PadPlane distance  */
  double innerSectorGatingGridPadSep; /*   GatingGrid-to-PadPlane distance  */
  double outerSectorAnodeWirePadSep; /*   AnodeWire-to-PadPlane distance  */
  double outerSectorFrischGridPadSep; /*   FrischGrid-to-PadPlane distance  */
  double outerSectorGatingGridPadSep; /*   GatingGrid-to-PadPlane distance  */
  int numInnerSectorAnodeWires; /*   */
  int numInnerSectorFrischGridWires; /*   */
  int numInnerSectorGatingGridWires; /*   */
  double firstInnerSectorAnodeWire; /*   */
  double firstInnerSectorFrischGridWire; /*   */
  double firstInnerSectorGatingGridWire; /*   */
  double lastInnerSectorAnodeWire; /*   */
  int numOuterSectorAnodeWires; /*   */
  int numOuterSectorFrischGridWires; /*   */
  int numOuterSectorGatingGridWires; /*   */
  double firstOuterSectorAnodeWire; /*   */
  double firstOuterSectorFrischGridWire; /*   */
  double firstOuterSectorGatingGridWire; /*   */
  double lastOuterSectorAnodeWire; /*   */
};

/** trigger offset common to all detector-level clock modules  trigger offset common to all detector-level clock modules */
struct trgTimeOffset {
  float offset; /*   standard trigger offset in micro-seconds  */
  float laserOffset; /*   laser trigger offset in micro-seconds  */
  float laserOffsetW; /*   laser extra trigger offset for West laser */
};

struct trigDetSums {
  unsigned int runNumber; /*       run number  */
  unsigned int timeOffset; /*       run begin time  */
  double ctbWest; /*   ctb West  */
  double ctbEast; /*   ctb East  */
  double ctbTOFp; /*   ctbOr + TOFp rate  */
  double tofp; /*   TOFp rate  */
  double zdcWest; /*    zdc west rate  */
  double zdcEast; /*    zdc east rate  */
  double zdcX; /*   zdc and rate  */
  double mult; /*   mult rate  */
  double L0; /*   L0 Rate  */
  double bbcX; /*   BBC and Rate  */
  double bbcXctbTOFp; /*   BBCAnd + ctbTOFp rate  */
  double bbcWest; /*   --BBC West--  */
  double bbcEast; /*   --BBC East--  */
  double bbcYellowBkg; /*   --(BBC Eastdelayed) and (BBC West)--  */
  double bbcBlueBkg; /*   --(BBC Westdelayed) and (BBC East)--  */
  double pvpdWest; /*   --PVPD East--  */
  double pvpdEast; /*   --PVPD West--  */
};

/** Table containg parameters for slow simulator running.  */
struct tss_tsspar {
  char fileout[80]; /* output file for pixel data (none->table) */
  int dynam; /* adc dynamic range (adc counts; usu 1023) */
  int format; /* pixel data format */
  int max_itime; /* upper bound of time bucket */
  int max_pads; /* upper bound of pads */
  int max_row; /* upper bound of row (<=45) */
  int max_sect; /* upper bound of sector pair (<=12) */
  int min_itime; /* low bound of time bucket */
  int min_pads; /* lower bound of pads */
  int min_row; /* lower bound of row (>=1) */
  int min_sect; /* lower bound of sector pair (>=1) */
  int mode; /* mode of TPC simulation for diff. tests */
  int nele_laser; /* number of electrons from laser point*/
  int ngain; /* number of gain sampling (1-10) */
  int nseg; /* number of sub-segment within a G-vol */
  int ntime; /* number of time buckets (dimensionless) */
  int printout; /* control the level of printout (0=no) */
  int tpc_half; /* half (1) or full (0) TPC volume */
  int reset; /* re-do setup: 0=no, anything else=yes */
  float ave_ion_pot; /* Average Ion. Potential of a gas(Ar=26eV) */
  float bfield; /* magnetic field strength (Tesla) */
  float c_test; /* test capacitance value (pF) */
  float diff_long; /* long diff const of gas (cm/sqrt(cm)) */
  float diff_trans; /* trans diff const of gas (cm/sqrt(cm)) */
  float gain_in; /* gas gain:inner sector (dimensionless) */
  float gain_out; /* gas gain:outer sector (dimensionless) */
  float prf_in; /* pad resp.func:inner sector (cm) */
  float prf_out; /* pad resp.func:outer sector (cm) */
  float sca_rms; /* SCA noise (random, not filtered) */
  float scale; /* number of electrons per ADC count (600) */
  float step_size; /* step size for subpadrow tracking */
  float tau; /* shaper resp. time const. (usec) */
  float threshold; /* adc threshold (adc counts but a real num */
  float time_offset; /* Wayne Betts time offsets */
  float v_test; /* voltage of test pulse (V) */
  float white_rms; /* rms of white noise (shaper filtered) */
  float wire_coupling_in; /* wire-to-pad coupling inner sector */
  float wire_coupling_out; /* wire-to-pad coupling outer sector */
  float x_laser; /* local x of laser point[cm] along row */
  float y_laser; /* local y of laser point[cm] across row */
  float z_laser; /* z drift length of pointlaser source(cm) */
};

#endif
