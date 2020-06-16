#ifndef TPCRS_STRUCT_CONTAINERS_H_
#define TPCRS_STRUCT_CONTAINERS_H_

#include "TArrayF.h"
#include "TArrayD.h"
#include "TGeoMatrix.h"
#include "TMath.h"
#include "TF1.h"

#include "tpcrs/configurator.h"
#include "tpcrs/structs.h"
#include "enums.h"
#include "config_structs.h"


struct St_itpcPadGainT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_itpcPadGainT0C, itpcPadGainT0>
{
  int 	run(int i = 0) 	const {return Struct(i)->run;}
  float 	Gain(int sector, int row, int pad) const
  {
    return ((sector > 0 && sector <= 24) && (row > 0 && row <= 40) && (pad > 0 && pad <= 120)) ?
           Struct()->Gain[sector - 1][row - 1][pad - 1] : 0;
  }
  float 	  T0(int sector, int row, int pad) const
  {
    return ((sector > 0 && sector <= 24) && (row > 0 && row <= 40) && (pad > 0 && pad <= 120)) ?
           Struct()->T0[sector - 1][row - 1][pad - 1] : 0;
  }
  bool    livePadrow(int sector, int row)
  {
    for (int pad = 1; pad <= 120; pad++) if (Gain(sector, row, pad) > 0) return kTRUE;

    return kFALSE;
  }
};

struct St_iTPCSurveyC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_iTPCSurveyC, iTPCSurvey> {
  int 	Id(int i = 0) 	const {return Struct(i)->Id;}
  float 	Angle(int i = 0) 	const {return Struct(i)->Angle;}
  float 	dx(int i = 0) 	const {return Struct(i)->dx;}
  float 	dy(int i = 0) 	const {return Struct(i)->dy;}
  float 	ScaleX(int i = 0) 	const {return Struct(i)->ScaleX;}
  float 	ScaleY(int i = 0) 	const {return Struct(i)->ScaleY;}
  char* 	comment(int i = 0) 	const {return Struct(i)->comment;}
};

struct St_MagFactorC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_MagFactorC, MagFactor> {
  float 	ScaleFactor(int i = 0) {return Struct(i)->ScaleFactor;}
};

struct St_MDFCorrectionC : tpcrs::IConfigStruct {
  virtual MDFCorrection* Struct(int i = 0) const = 0;
  enum EMDFPolyType {
    kMonomials,
    kChebyshev,
    kLegendre
  };
  St_MDFCorrectionC();

  void Initialize()
  {
    unsigned int N = GetNRows();
    fFunc = new TF1*[N];
    memset(fFunc, 0, N*sizeof(TF1*));
  }

  unsigned char 	idx(int k = 0)        	const {return Struct(k)->idx;}
  unsigned char 	nrows(int k = 0) 	        const {return Struct(k)->nrows;}
  unsigned char 	PolyType(int k = 0) 	        const {return Struct(k)->PolyType;}
  unsigned char 	NVariables(int k = 0) 	const {return Struct(k)->NVariables;}
  unsigned char 	NCoefficients(int k = 0) 	const {return Struct(k)->NCoefficients;}
  unsigned char* 	Powers(int k = 0) 	        const {return Struct(k)->Power;}
  double 	DMean(int k = 0)           	const {return Struct(k)->DMean;}
  double* 	XMin(int k = 0)         	const {return Struct(k)->XMin;}
  double* 	XMax(int k = 0)       	const {return Struct(k)->XMax;}
  double* 	Coefficients(int k = 0) 	const {return Struct(k)->Coefficients;}
  double* 	CoefficientsRMS(int k = 0) 	const {return Struct(k)->CoefficientsRMS;}
  double      Eval(int k = 0, double *x = 0) const;
  double      Eval(int k, double x0, double x1) const;
  double      EvalError(int k = 0, double *x = 0) const;
  static double MDFunc(double *x = 0, double *p = 0);
  static St_MDFCorrectionC *fgMDFCorrectionC;
 protected:
  virtual ~St_MDFCorrectionC();
 private:
  double EvalFactor(int k = 0, int p = 0, double x = 0) const;
  TF1         **fFunc;
};

struct St_richvoltagesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_richvoltagesC, richvoltages>
{
  unsigned int 	runNumber(int i = 0) 	        {return Struct(i)->runNumber;}
  unsigned int 	startStatusTime(int i = 0) 	{return Struct(i)->startStatusTime;}
  unsigned int 	endStatusTime(int i = 0) 	{return Struct(i)->endStatusTime;}
  unsigned int 	status(int i = 0) 	        {return Struct(i)->status;}
};

enum StMagnetPolarity {eUnknownMField, eFullMFieldPolB, eHalfMFieldPolB,
                       eZeroMField, eHalfMFieldPolA, eFullMFieldPolA
                      };

struct St_starMagOnlC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_starMagOnlC, starMagOnl>
{
  unsigned int 	runNumber(int i = 0) 	{return Struct(i)->runNumber;}
  unsigned int 	time(int i = 0) 	{return Struct(i)->time;}
  double 	current(int i = 0) 	{return Struct(i)->current;}
  double      getScaleFactor(unsigned int time = 0) {return currentToScaleFactor(getMagnetCurrent(time));}
  double      getMagnetCurrent(unsigned int time = 0)
  {
    if (GetNRows() == 1 || time == 0) return current();

    double tempCurrent = -9999;

    for (unsigned int i = 0; i < GetNRows() - 1; i++)
      if ( time >= getTimeEntry(i) && time <= getTimeEntry(i + 1) )
        if ( TMath::Abs(getMagnetCurrentEntry(i) - getMagnetCurrentEntry(i + 1)) < 50 )
          tempCurrent = getMagnetCurrentEntry(i);

    return tempCurrent;
  }
  StMagnetPolarity           getMagneticField(unsigned int time = 0)
  {
    StMagnetPolarity value = eUnknownMField;

    double scaleFactor = getScaleFactor(time);

    if (scaleFactor == 1.0)	value = eFullMFieldPolA;

    if (scaleFactor == 0.5)	value = eHalfMFieldPolA;

    if (scaleFactor == 0.0)	value = eZeroMField;

    if (scaleFactor == -0.5)	value = eHalfMFieldPolB;

    if (scaleFactor == -1.0)	value = eFullMFieldPolB;

    return value;
  }
  unsigned int        getRunNumber() {return runNumber();}
  unsigned int        getTimeEntry(unsigned int i = 0) {return time(i);}
  double      getMagnetCurrentEntry(unsigned int i = 0) {return current(i);}
  static double  currentToScaleFactor(double current)
  {
    double value = -9999;

    if     (current < -4450 && current > -4550)	value = -1.0;
    else if (current < -2200 && current > -2300)	value = -0.5;
    else if (current >   -50 && current <    50)	value =  0.0;
    else if (current >  2200 && current <  2300)	value =  0.5;
    else if (current >  4450 && current <  4550)	value =  1.0;

    return value;
  }
};

struct St_spaceChargeCorC : tpcrs::IConfigStruct {
  virtual spaceChargeCor* Struct(int i = 0) const = 0;
  double 	fullFieldB(int i = 0) 	{return Struct(i)->fullFieldB;}
  double 	halfFieldB(int i = 0) 	{return Struct(i)->halfFieldB;}
  double 	zeroField(int i = 0) 	        {return Struct(i)->zeroField;}
  double 	halfFieldA(int i = 0) 	{return Struct(i)->halfFieldA;}
  double 	fullFieldA(int i = 0) 	{return Struct(i)->fullFieldA;}
  double 	satRate(int i = 0) 	        {return Struct(i)->satRate;}
  float 	factor(int i = 0) 	        {return Struct(i)->factor;}
  float 	detector(int i = 0) 	        {return Struct(i)->detector;}
  float 	offset(int i = 0) 	        {return Struct(i)->offset;}
  float 	getEWRatio(int i = 0)	        {return Struct(i)->ewratio;}
  double      getSpaceChargeCorrection(double scaleFactor, int i = 0){
    double value = 0;
    if(scaleFactor < -.75 && scaleFactor > -1.25) value = fullFieldB(i);
    else if(scaleFactor < -0.25)	          value = halfFieldB(i);
    else if(scaleFactor < .25)	                  value = zeroField(i);
    else if(scaleFactor < 0.75)	                  value = halfFieldA(i);
    else if(scaleFactor < 1.25)	                  value = fullFieldA(i);
    return value;
  }
  double getSpaceChargeCorrection(){return  getSpaceChargeCorrection(St_starMagOnlC::instance()->getScaleFactor());}
  double getSpaceChargeCoulombs(double scaleFactor);
  double getSpaceChargeCoulombs(){return getSpaceChargeCoulombs(St_starMagOnlC::instance()->getScaleFactor());}
  double getSpaceChargeSatRate(int i = 0) {return satRate(i);}
  float  getSpaceChargeFactor(int i = 0)  {return factor(i);}
  float  getSpaceChargeDetector(int i = 0){return detector(i);}
  float  getSpaceChargeOffset(int i = 0)  {return offset(i);}

};

struct St_spaceChargeCorR1C : tpcrs::ConfigStruct<St_spaceChargeCorC, St_spaceChargeCorR1C, spaceChargeCor> {};

struct St_spaceChargeCorR2C : tpcrs::ConfigStruct<St_spaceChargeCorC, St_spaceChargeCorR2C, spaceChargeCor> {};

struct St_SurveyC : tpcrs::IConfigStruct {
  virtual Survey* Struct(int i = 0) const = 0;

  void Initialize()
  {
    unsigned int N = GetNRows();
    fRotations = new TGeoHMatrix*[N];
    for (unsigned int i = 0; i < N; i++) {
      fRotations[i] = new TGeoHMatrix;
      TGeoHMatrix &rot = *fRotations[i];
      if (N == 1) rot.SetName(GetName().c_str());
      else        rot.SetName(Form("%s_%i",GetName().c_str(),i+1));
      rot.SetRotation(Rotation(i));
      rot.SetTranslation(Translation(i));
      Normalize(rot);
      assert(TMath::Abs(rot.Determinant())-1 < 1.e-3);
    }
  }

  virtual  ~St_SurveyC();
  int 	Id(int i = 0) 	const {return Struct(i)->Id;}
  double 	r00(int i = 0) 	const {return Struct(i)->r00;} // 0
  double 	r01(int i = 0) 	const {return Struct(i)->r01;} // 1
  double 	r02(int i = 0) 	const {return Struct(i)->r02;} // 2
  double 	r10(int i = 0) 	const {return Struct(i)->r10;} // 3
  double 	r11(int i = 0) 	const {return Struct(i)->r11;} // 4
  double 	r12(int i = 0) 	const {return Struct(i)->r12;} // 5
  double 	r20(int i = 0) 	const {return Struct(i)->r20;} // 6
  double 	r21(int i = 0) 	const {return Struct(i)->r21;} // 7
  double 	r22(int i = 0) 	const {return Struct(i)->r22;} // 8
  double 	t0(int i = 0) 	const {return Struct(i)->t0;}
  double 	t1(int i = 0) 	const {return Struct(i)->t1;}
  double 	t2(int i = 0) 	const {return Struct(i)->t2;}
  double 	sigmaRotX(int i = 0) 	const {return Struct(i)->sigmaRotX;}
  double 	sigmaRotY(int i = 0) 	const {return Struct(i)->sigmaRotY;}
  double 	sigmaRotZ(int i = 0) 	const {return Struct(i)->sigmaRotZ;}
  double 	sigmaTrX(int i = 0) 	const {return Struct(i)->sigmaTrX;}
  double 	sigmaTrY(int i = 0) 	const {return Struct(i)->sigmaTrY;}
  double 	sigmaTrZ(int i = 0) 	const {return Struct(i)->sigmaTrZ;}
  char* 	comment(int i = 0) 	const {return Struct(i)->comment;}
  void          GetAngles(double &phi, double &the, double &psi, int i = 0);
  const double  *Rotation(int i = 0)     const {return &Struct(i)->r00;} 
  const double  *Translation(int i = 0)  const {return &Struct(i)->t0;} 
  const TGeoHMatrix  &GetMatrix(int i = 0);
  const TGeoHMatrix  &GetMatrix4Id(int id);
  const TGeoHMatrix  &GetMatrixR(int i); // ignoring rotation alpha and beta
  const double *r(int i = 0)        const {return &Struct(i)->r00;}
  const double *t(int i = 0)        const {return &Struct(i)->t0;}
  static void Normalize(TGeoHMatrix &rot);
  static double IsOrtogonal(const double *r);
 protected:
  St_SurveyC();
 private:
  TGeoHMatrix  **fRotations;
};

struct St_TpcAdcCorrectionMDF : tpcrs::ConfigStruct<St_MDFCorrectionC, St_TpcAdcCorrectionMDF, MDFCorrection> {};

struct St_tpcAnodeHVavgC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAnodeHVavgC, tpcAnodeHVavg>
{
  unsigned short          sector(int i = 0) 	const {return Struct(i)->sector;}
  unsigned short          socket(int i = 0) 	const {return Struct(i)->socket;}
  float 	    voltage(int i = 0) 	const;
  float 	    rms(int i = 0) 	        const {return Struct(i)->rms;}
  int 	    numentries(int i = 0) 	const {return Struct(i)->numentries;}
  int 	    numoutliers(int i = 0) 	const {return Struct(i)->numoutliers;}
  bool	    livePadrow(int sec = 1, int padrow = 1) const { return voltagePadrow(sec, padrow) > 500; }
  float	    voltagePadrow(int sec = 1, int padrow = 1) const; // sector=1..24 , padrow=1..100
  bool            tripped(int sec = 1, int padrow = 1)       const;// { return (voltage() < -100); }
};

struct St_tpcAnodeHVC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcAnodeHVC, tpcAnodeHV>
{
  unsigned short 	 sector(int i = 0) 	const {return Struct(i)->sector;}
  unsigned short 	 socket(int i = 0) 	const {return Struct(i)->socket;}
  float 	 voltage(int i = 0) 	const;
  bool	 livePadrow(int sector = 1, int padrow = 1) const { return voltagePadrow(sector, padrow) > 500; }
  float	 voltagePadrow(int sector = 1, int padrow = 1) const ; // sector=1..24 , padrow=1..100
  bool         tripped(int sector = 1, int padrow = 1) const { return (voltagePadrow(sector, padrow) < -100); }
  static  void   sockets(int sector, int padrow, int &e1, int &e2, float &f2);
};

struct St_TpcAvgCurrentC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcAvgCurrentC, TpcAvgCurrent>
{
  int 	run(int i = 0) 	const {return Struct(i)->run;}
  int 	start_time(int i = 0) 	const {return Struct(i)->start_time;}
  int 	stop_time(int i = 0) 	const {return Struct(i)->stop_time;}
  static int  ChannelFromRow(int sector, int row);
  static int  ChannelFromSocket(int socket);
  float       AvCurrent(int sector = 1, int channel = 1);
  /* {
     return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
     Struct()->AvCurrent[8*(sector-1)+channel-1] :
     0;} */
  float       AvCurrSocket(int sector = 1, int socket = 1) {return AvCurrent(sector, ChannelFromSocket(socket));}
  float       AvCurrRow(int sector = 1, int row = 1) {return AvCurrent(sector, ChannelFromRow(sector, row));}
  float       AcCharge(int sector = 1, int channel = 1);
  /* {
     return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
     Struct()->AcCharge[8*(sector-1)+channel-1] :
     0;
     } */
  float       AcChargeSocket(int sector = 1, int socket = 1) {return AcCharge(sector, ChannelFromSocket(socket));}
  float       AcChargeRow(int sector = 1, int row = 1) {return AcCharge(sector, ChannelFromRow(sector, row));}
  float       AcChargeL(int sector = 1, int channel = 1); // C/cm
  float       AcChargeRowL(int sector = 1, int row = 1) {return AcChargeL(sector, ChannelFromRow(sector, row));}
};

struct St_TpcAvgPowerSupplyC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcAvgPowerSupplyC, TpcAvgPowerSupply>
{
  int 	run(int i = 0) 	const {return Struct(i)->run;}
  int 	start_time(int i = 0) 	const {return Struct(i)->start_time;}
  int 	stop_time(int i = 0) 	const {return Struct(i)->stop_time;}
  float* 	Current(int i = 0) 	const {return Struct(i)->Current;}
  float* 	Charge(int i = 0) 	const {return Struct(i)->Charge;}
  float* 	Voltage(int i = 0) 	const {return Struct(i)->Voltage;}
  float	voltagePadrow(int sec = 1, int padrow = 1) const; // sector=1..24 , padrow=1..100
  bool        tripped(int sec = 1, int row = 1) const {return voltagePadrow(sec, row) < -100;}
  static int  ChannelFromRow(int sector, int row) {return St_TpcAvgCurrentC::ChannelFromRow(sector, row);}
  static int  ChannelFromSocket(int socket) {return St_TpcAvgCurrentC::ChannelFromSocket(socket);}
  float       AvCurrent(int sector = 1, int channel = 1)
  {
    return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
           Struct()->Current[8 * (sector - 1) + channel - 1] :
           0;
  }
  float       AvCurrSocket(int sector = 1, int socket = 1) {return AvCurrent(sector, ChannelFromSocket(socket));}
  float       AvCurrRow(int sector = 1, int row = 1) {return AvCurrent(sector, ChannelFromRow(sector, row));}
  float       AcCharge(int sector = 1, int channel = 1)
  {
    return (sector > 0 && sector <= 24 && channel > 0 && channel <= 8) ?
           Struct()->Charge[8 * (sector - 1) + channel - 1] :
           0;
  }
  float       AcChargeSocket(int sector = 1, int socket = 1) {return AcCharge(sector, ChannelFromSocket(socket));}
  float       AcChargeRow(int sector = 1, int row = 1) {return AcCharge(sector, ChannelFromRow(sector, row));}
  float       AcChargeL(int sector = 1, int channel = 1); // C/cm
  float       AcChargeRowL(int sector = 1, int row = 1) {return AcChargeL(sector, ChannelFromRow(sector, row));}
  bool        livePadrow(int sec = 1, int padrow = 1) const { return voltagePadrow(sec, padrow) >  500;}
};

struct St_tpcCalibResolutionsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcCalibResolutionsC, tpcCalibResolutions> {
  float 	SpaceCharge(int i = 0) 	{return Struct(i)->SpaceCharge;}
  float 	GridLeak(int i = 0)	 	{return Struct(i)->GridLeak;}
};

struct St_tpcChargeEventC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcChargeEventC, tpcChargeEvent> {
  int nChargeEvents()                            {return Struct()->nChargeEvents;}
  unsigned int* eventBunchCrossingsLow()         {return Struct()->eventBunchCrossingsLow;}
  unsigned int* eventBunchCrossingsHigh()        {return Struct()->eventBunchCrossingsHigh;}
  float* eventCharges()                          {return Struct()->eventCharges;}
  int badBunch()                                 {return Struct()->badBunch;}

  unsigned int eventBunchCrossingLow(int idx)    {return eventBunchCrossingsLow()[idx]; }
  unsigned int eventBunchCrossingHigh(int idx)   {return eventBunchCrossingsHigh()[idx]; }
  unsigned long long eventBunchCrossing(int idx) {return (((unsigned long long) (eventBunchCrossingHigh(idx))) << 32)
                                                        + ((unsigned long long) (eventBunchCrossingLow(idx))); }
  float eventCharge(int idx)                     {return eventCharges()[idx];}

  // user functions for getting the charge and time since charge

  void lastChargeTime(unsigned long long bunchCrossingNumber, float& charge, double& timeSinceCharge) {
    int idx = indexBeforeBunchCrossing(bunchCrossingNumber);
    charge = eventCharge(idx);
    timeSinceCharge = timeDifference(bunchCrossingNumber,idx);
  }

  // must call findLastChargeTime() before getLastChargeTime()
  void findLastchargeTime(unsigned long long bunchCrossingNumber) {
    lastChargeTime(bunchCrossingNumber, localStoreCharge, localStoreTimeSinceCharge);
  }
  void getLastChargeTime(float& charge, double& timeSinceCharge) {
    charge = localStoreCharge;
    timeSinceCharge = localStoreTimeSinceCharge;
  }

  // must call findChargeTimes() before getCharges() and getTimes()
  int findChargeTimes(unsigned long long bunchCrossingNumber, unsigned long long bunchCrossingWindow);
  int findChargeTimes(unsigned long long bunchCrossingNumber, double timeWindow=1.9);
  TArrayF* getCharges() { return &localStoreCharges; }
  TArrayD* getTimes() { return &localStoreTimesSinceCharges; }

 protected:
  double timeDifference(unsigned long long bunchCrossingNumber, int idx);
  int indexBeforeBunchCrossing(unsigned long long bunchCrossingNumber);
 private:
  int localSearchLowerIndex = 0;
  int localSearchUpperIndex = -1;
  float localStoreCharge = 0; //!
  double localStoreTimeSinceCharge = 0; //!
  TArrayF localStoreCharges; //!
  TArrayD localStoreTimesSinceCharges; //!
};

struct St_tpcCorrectionC : tpcrs::IConfigStruct
{
  virtual tpcCorrection* Struct(int i = 0) const = 0;
  int 	type(int i = 0) 	const {return Struct(i)->type;}
  int 	idx(int i = 0) 	const {return Struct(i)->idx;}
  int 	nrows(int i = 0) 	const {return Struct(i)->nrows;}
  int 	npar(int i = 0) 	const {return Struct(i)->npar;}
  double 	OffSet(int i = 0) 	const {return Struct(i)->OffSet;}
  double 	min(int i = 0) 	const {return Struct(i)->min;}
  double 	max(int i = 0) 	const {return Struct(i)->max;}
  double* 	a(int i = 0) 	        const {return Struct(i)->a;}
  double CalcCorrection(int i, double x, double z = 0, int NparMax = -1);
  double SumSeries(tpcCorrection* cor, double x, double z = 0, int NparMax = -1);
};

struct St_TpcAdcCorrectionBC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcAdcCorrectionBC, tpcCorrection> {};

struct St_TpcCurrentCorrectionC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcCurrentCorrectionC, tpcCorrection> {};

struct St_TpcdChargeC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcdChargeC, tpcCorrection> {};

struct St_TpcdEdxCorC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcdEdxCorC, tpcCorrection> {};

struct St_tpcPadConfigC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadConfigC, tpcPadConfig>
{
  unsigned char          iTpc(int sector);
  unsigned char          iTPC(int sector) {return iTpc(sector);}
  int 	   padRows(int sector);
  int 	   innerPadRows(int sector);
  int 	   innerPadRows48(int sector);
  int 	   innerPadRows52(int sector);
  int 	   outerPadRows(int sector);
  int 	   superInnerPadRows(int sector);
  int 	   superOuterPadRows(int sector);
  double 	   innerSectorPadWidth(int sector);
  double 	   innerSectorPadLength(int sector);
  double 	   innerSectorPadPitch(int sector);
  double 	   innerSectorRowPitch1(int sector);
  double 	   innerSectorRowPitch2(int sector);
  double 	   firstPadRow(int sector);
  double 	   firstOuterSectorPadRow(int sector);
  double 	   lastOuterSectorPadRow(int sector);
  double 	   firstRowWidth(int sector);
  double 	   lastRowWidth(int sector);
  double 	   outerSectorPadWidth(int sector);
  double 	   outerSectorPadLength(int sector);
  double 	   outerSectorPadPitch(int sector);
  double 	   outerSectorRowPitch(int sector);
  double 	   outerSectorLength(int sector);
  double 	   ioSectorSeparation(int sector);
  double 	   innerSectorEdge(int sector);
  double 	   outerSectorEdge(int sector);
  double 	   innerSectorPadPlaneZ(int sector);
  double 	   outerSectorPadPlaneZ(int sector);
  int* 	   innerPadsPerRow(int sector);
  int* 	   outerPadsPerRow(int sector);
  int            padsPerRow(int sector, int row = 1);
  double* 	   innerRowRadii(int sector);
  double* 	   outerRowRadii(int sector);
  //               taken from StRItpcPadPlane
  int            numberOfRows(int sector);
  int            numberOfInnerRows(int sector);
  int            numberOfInnerRows48(int sector);
  int            numberOfInnerRows52(int sector);
  int            numberOfOuterRows(int sector);
  bool           isRowInRange(int sector, int row);
  double         radialDistanceAtRow(int sector, int row);
  int            numberOfPadsAtRow(int sector, int row);
  double         PadWidthAtRow(int sector, int row);
  double 	   PadLengthAtRow(int sector, int row);
  double 	   PadPitchAtRow(int sector, int row);
  double 	   RowPitchAtRow(int sector, int row);
  int            indexForRowPad(int sector, int row, int pad);
  bool             isiTpcSector(int sector) { return iTpc(sector) == 1; }
  bool             isiTpcPadRow(int sector, int row) { return iTpc(sector) && row >= 1 && row <= numberOfInnerRows(sector); }
  bool             isInnerPadRow(int sector, int row) { return row <= numberOfInnerRows(sector); }
  int            IsRowInner(int sector, int row) {return (row <= innerPadRows(sector)) ? 1 : 0;}
};

struct St_TpcDriftDistOxygenC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcDriftDistOxygenC, tpcCorrection> {};

struct St_tpcDriftVelocityC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcDriftVelocityC, tpcDriftVelocity> {
  float 	laserDriftVelocityEast(int i = 0) 	{return Struct(i)->laserDriftVelocityEast;}
  float 	laserDriftVelocityWest(int i = 0) 	{return Struct(i)->laserDriftVelocityWest;}
  float 	cathodeDriftVelocityEast(int i = 0) 	{return Struct(i)->cathodeDriftVelocityEast;}
  float 	cathodeDriftVelocityWest(int i = 0) 	{return Struct(i)->cathodeDriftVelocityWest;}
};

struct St_TpcdXCorrectionBC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcdXCorrectionBC, tpcCorrection> {};

struct St_TpcEdgeC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcEdgeC, tpcCorrection> {};

struct St_TpcEffectivedXC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_TpcEffectivedXC, TpcEffectivedX> {
  float 	scaleInner(int i = 0) 	const {return Struct(i)->scaleInner;}
  float 	scaleOuter(int i = 0) 	const {return Struct(i)->scaleOuter;}
};

struct St_tpcFieldCageC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcFieldCageC, tpcFieldCage>
{
  float 	innerFieldCageShift(int i = 0) {return Struct(i)->innerFieldCageShift;}
  float 	InnerFieldCageShift(int i = 0) {return innerFieldCageShift(i);}
  float 	eastClockError(int i = 0) 	{return Struct(i)->eastClockError;}
  float 	EastClockError(int i = 0) 	{return eastClockError(i);}
  float 	westClockError(int i = 0) 	{return Struct(i)->westClockError;}
  float 	WestClockError(int i = 0) 	{return westClockError(i);}
};

struct St_tpcFieldCageShortC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcFieldCageShortC, tpcFieldCageShort> {
  float 	side(int i = 0) 	        {return Struct(i)->side;}
  float 	cage(int i = 0) 	        {return Struct(i)->cage;}
  float 	ring(int i = 0) 	        {return Struct(i)->ring;}
  float 	resistor(int i = 0) 	        {return Struct(i)->resistor;}
  float 	MissingResistance(int i = 0) 	{return Struct(i)->MissingResistance;}
};

struct St_tpcGainCorrectionC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_tpcGainCorrectionC, tpcCorrection>
{};

struct St_tpcGasTemperatureC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_tpcGasTemperatureC, tpcCorrection> {};

struct St_tpcGlobalPositionC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGlobalPositionC, tpcGlobalPosition>
{
  float 	LocalxShift(int i = 0)       const {return Struct(i)->LocalxShift;}
  float 	LocalyShift(int i = 0)       const {return Struct(i)->LocalyShift;}
  float 	LocalzShift(int i = 0)       const {return Struct(i)->LocalzShift;}
  /*  float 	PhiXY(int i = 0)  	       const {return Struct(i)->PhiXY;}	   */
  float 	PhiXZ(int i = 0)  	       const {return Struct(i)->PhiXZ;}
  float 	PhiYZ(int i = 0)  	       const {return Struct(i)->PhiYZ;}
  /*  float 	XX(int i = 0)  	       const {return Struct(i)->XX;}
      float 	YY(int i = 0)  	       const {return Struct(i)->YY;}
      float 	ZZ(int i = 0)  	       const {return Struct(i)->ZZ;}	    */
  float 	PhiXY_geom(int i = 0)        const {return Struct(i)->PhiXY_geom;}
  float 	PhiXZ_geom(int i = 0)        const {return Struct(i)->PhiXZ_geom;}
  float 	PhiYZ_geom(int i = 0)        const {return Struct(i)->PhiYZ_geom;}
  /*  float 	XX_geom(int i = 0)  	       const {return Struct(i)->XX_geom;}
      float 	YY_geom(int i = 0)  	       const {return Struct(i)->YY_geom;}
      float 	ZZ_geom(int i = 0)  	       const {return Struct(i)->ZZ_geom;}   */
  double  	TpcCenterPositionX()           const {return LocalxShift();}
  double  	TpcCenterPositionY()           const {return LocalyShift();}
  double  	TpcCenterPositionZ()           const {return LocalzShift();}
  double  	TpcRotationAroundGlobalAxisX() const {return PhiYZ_geom();}
  double  	TpcRotationAroundGlobalAxisY() const {return PhiXZ_geom();}
  double  	TpcRotationAroundGlobalAxisZ() const {return PhiXY_geom();}
  double  	TpcEFieldRotationX()           const {return PhiYZ();} /* YTWIST */
  double  	TpcEFieldRotationY() 	       const {return PhiXZ();} /* XTWIST */
  double      XTWIST()                       const {return  1e3 * TpcEFieldRotationY();}
  double      YTWIST()                       const {return -1e3 * TpcEFieldRotationX();}
  /* double  	TpcEFieldRotationZ() 	       const {return PhiXY();}              */
  double      X0()                           const {return LocalxShift();}
  double      Y0()                           const {return LocalyShift();}
  double      Z0()                           const {return LocalzShift();}
  double      alpha()                        const {return PhiYZ_geom();}
  double      beta()                         const {return PhiXZ_geom();}
  double      gamma()                        const {return PhiXY_geom();}
};

enum StGLpos {
  kGLinner=0,
  kGLmiddl=1,
  kGLouter=2
};

struct St_tpcGridLeakC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcGridLeakC, tpcGridLeak> {
  double 	InnerGLRadius(int i = 0) 	{return Struct(i)->InnerGLRadius;}
  double 	MiddlGLRadius(int i = 0) 	{return Struct(i)->MiddlGLRadius;}
  double 	OuterGLRadius(int i = 0) 	{return Struct(i)->OuterGLRadius;}
  double 	InnerGLWidth(int i = 0) 	{return Struct(i)->InnerGLWidth;}
  double 	MiddlGLWidth(int i = 0) 	{return Struct(i)->MiddlGLWidth;}
  double 	OuterGLWidth(int i = 0) 	{return Struct(i)->OuterGLWidth;}
  double 	InnerGLStrength(int i = 0) 	{return Struct(i)->InnerGLStrength;}
  double 	MiddlGLStrength(int i = 0) 	{return Struct(i)->MiddlGLStrength;}
  double 	OuterGLStrength(int i = 0) 	{return Struct(i)->OuterGLStrength;}
  double      getGridLeakStrength(StGLpos pos){
    switch (pos) {
      case (kGLinner) : return InnerGLStrength();
      case (kGLmiddl) : return MiddlGLStrength();
      case (kGLouter) : return OuterGLStrength();
    }
    return 0;
  }
  double      getGridLeakRadius(StGLpos pos) {
    switch (pos) {
      case (kGLinner) : return InnerGLRadius();
      case (kGLmiddl) : return MiddlGLRadius();
      case (kGLouter) : return OuterGLRadius();
    }
    return 0;
  }
  double      getGridLeakWidth(StGLpos pos) {
    switch (pos) {
      case (kGLinner) : return InnerGLWidth();
      case (kGLmiddl) : return MiddlGLWidth();
      case (kGLouter) : return OuterGLWidth();
    }
    return 0;
  }
};

struct St_tpcHighVoltagesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcHighVoltagesC, tpcHighVoltages> {
  float 	cathode(int i = 0)          {return Struct(i)->cathode;}
  float 	gatedGridRef(int i = 0)     {return Struct(i)->gatedGridRef;}
  float* 	gridLeakWallTip(int i = 0)  {return Struct(i)->gridLeakWallTip;}
  float* 	gridLeakWallSide(int i = 0) {return Struct(i)->gridLeakWallSide;}
  double      getCathodeVoltage()           {return cathode();}
  double      getGGVoltage()                {return gatedGridRef();}
  double      getGridLeakWallTip(int sector = 1)  {return gridLeakWallTip()[sector-1];}
  double      getGridLeakWallSide(int sector = 1) {return gridLeakWallSide()[sector-1];}
};

struct St_tpcHVPlanesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcHVPlanesC, tpcHVPlanes> {};

struct St_TpcLengthCorrectionBC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcLengthCorrectionBC, tpcCorrection> {};

struct St_TpcLengthCorrectionMDF : tpcrs::ConfigStruct<St_MDFCorrectionC, St_TpcLengthCorrectionMDF, MDFCorrection> {};

struct St_tpcMethaneInC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_tpcMethaneInC, tpcCorrection> {};

struct St_TpcMultiplicityC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcMultiplicityC, tpcCorrection> {};

struct St_tpcOmegaTauC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcOmegaTauC, tpcOmegaTau> {
 public:
  float 	tensorV1(int i = 0) 	     {return Struct(i)->tensorV1;}
  float 	tensorV2(int i = 0) 	     {return Struct(i)->tensorV2;}
  float 	getOmegaTauTensorV1()        {return tensorV1();}
  float 	getOmegaTauTensorV2()        {return tensorV2();}
  unsigned int        distortionCorrectionsMode(int i = 0)
                                             {return Struct(i)->distortionCorrectionsMode;}
};

struct St_TpcPadCorrectionMDF : tpcrs::ConfigStruct<St_MDFCorrectionC, St_TpcPadCorrectionMDF, MDFCorrection> {};

struct St_tpcPadGainT0BC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadGainT0BC, tpcPadGainT0>
{
  float 	Gain(int sector, int row, int pad) const;
  float 	  T0(int sector, int row, int pad) const;
  bool    livePadrow(int sector, int row) const;
};

struct St_tpcPadGainT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadGainT0C, tpcPadGainT0>
{
  int 	run()           	const {return Struct()->run;}
  float 	Gain(int sector, int row, int pad) const
  {
    float gain = 0;

    if ((sector > 0 && sector <= 24) && (row > 0 && row <= St_tpcPadConfigC::instance()->padRows(sector)) && (pad > 0 && pad <= 182)) {
      gain = Struct()->Gain[sector - 1][row - 1][pad - 1];
    }

    return gain;
  }
  float 	  T0(int sector, int row, int pad) const
  {
    float t0 = 0;

    if ((sector > 0 && sector <= 24) && (row > 0 && row <= St_tpcPadConfigC::instance()->padRows(sector)) && (pad > 0 && pad <= 182)) {
      t0 = Struct()->T0[sector - 1][row - 1][pad - 1];
    }

    return t0;
  }
  bool    livePadrow(int sector, int row)
  {
    for (int pad = 1; pad <= 182; pad++) if (Gain(sector, row, pad) > 0) return kTRUE;

    return kFALSE;
  }
};

struct St_tpcPadPlanesC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadPlanesC, tpcPadPlanes>
{
  int 	padRows(int i = 0) 	         {return Struct(i)->padRows;}
  int 	innerPadRows(int i = 0) 	 {return Struct(i)->innerPadRows;}
  int 	innerPadRows48(int i = 0) 	 {return Struct(i)->innerPadRows48;}
  int 	innerPadRows52(int i = 0) 	 {return Struct(i)->innerPadRows52;}
  int 	outerPadRows(int i = 0) 	 {return Struct(i)->outerPadRows;}
  int 	superInnerPadRows(int i = 0) 	 {return Struct(i)->superInnerPadRows;}
  int 	superOuterPadRows(int i = 0) 	 {return Struct(i)->superOuterPadRows;}
  double 	innerSectorPadWidth(int i = 0) {return Struct(i)->innerSectorPadWidth;}
  double 	innerSectorPadLength(int i = 0) {return Struct(i)->innerSectorPadLength;}
  double 	innerSectorPadPitch(int i = 0) {return Struct(i)->innerSectorPadPitch;}
  double 	innerSectorRowPitch1(int i = 0) {return Struct(i)->innerSectorRowPitch1;}
  double 	innerSectorRowPitch2(int i = 0) {return Struct(i)->innerSectorRowPitch2;}
  double 	firstPadRow(int i = 0) 	 {return Struct(i)->firstPadRow;}
  double 	firstOuterSectorPadRow(int i = 0) {return Struct(i)->firstOuterSectorPadRow;}
  double 	lastOuterSectorPadRow(int i = 0) {return Struct(i)->lastOuterSectorPadRow;}
  double 	firstRowWidth(int i = 0) 	 {return Struct(i)->firstRowWidth;}
  double 	lastRowWidth(int i = 0) 	 {return Struct(i)->lastRowWidth;}
  double 	outerSectorPadWidth(int i = 0) {return Struct(i)->outerSectorPadWidth;}
  double 	outerSectorPadLength(int i = 0) {return Struct(i)->outerSectorPadLength;}
  double 	outerSectorPadPitch(int i = 0) {return Struct(i)->outerSectorPadPitch;}
  double 	outerSectorRowPitch(int i = 0) {return Struct(i)->outerSectorRowPitch;}
  double 	outerSectorLength(int i = 0) 	 {return Struct(i)->outerSectorLength;}
  double 	ioSectorSeparation(int i = 0)  {return Struct(i)->ioSectorSeparation;}
  double 	innerSectorEdge(int i = 0) 	 {return Struct(i)->innerSectorEdge;}
  double 	outerSectorEdge(int i = 0) 	 {return Struct(i)->outerSectorEdge;}
  double 	innerSectorPadPlaneZ(int i = 0) {return Struct(i)->innerSectorPadPlaneZ;}
  double 	outerSectorPadPlaneZ(int i = 0) {return Struct(i)->outerSectorPadPlaneZ;}
  int* 	innerPadsPerRow(int i = 0) 	 {return Struct(i)->innerPadsPerRow;}
  int* 	outerPadsPerRow(int i = 0) 	 {return Struct(i)->outerPadsPerRow;}
  int         padsPerRow(int row = 1)
  {
    return (row <= innerPadRows()) ?
           innerPadsPerRow()[row - 1] :
           outerPadsPerRow()[row - 1 - innerPadRows()];
  }
  double* 	innerRowRadii(int i = 0) 	 {return Struct(i)->innerRowRadii;}
  double* 	outerRowRadii(int i = 0) 	 {return Struct(i)->outerRowRadii;}
  // taken from StRTpcPadPlane
  int         numberOfRows()                   {return padRows();}
  int         numberOfInnerRows()              {return innerPadRows();}
  int         numberOfInnerRows48()            {return innerPadRows48();}
  int         numberOfInnerRows52()            {return innerPadRows52();}
  int         numberOfOuterRows()              {return outerPadRows();}
  bool        isRowInRange(int row)          {return (row >= 1 && row <= numberOfRows()) ? kTRUE : kFALSE;}
  double      radialDistanceAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows() ) return innerRowRadii()[row - 1];
    else                            return outerRowRadii()[row - 1 - numberOfInnerRows()];
  }
  int   numberOfPadsAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows() ) return innerPadsPerRow()[row - 1];

    return outerPadsPerRow()[row - 1 - numberOfInnerRows()];
  }
  double PadWidthAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows()) return innerSectorPadWidth();

    return outerSectorPadWidth();
  }
  double PadLengthAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows()) return innerSectorPadLength();

    return outerSectorPadLength();
  }
  double PadPitchAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows()) return innerSectorPadPitch();

    return outerSectorPadPitch();
  }
  double RowPitchAtRow(int row)
  {
    if (! isRowInRange(row)) return 0;

    if ( row <= numberOfInnerRows48() ) return innerSectorRowPitch1();
    else if (row > numberOfInnerRows48() && row <= numberOfInnerRows()) return innerSectorRowPitch2();

    return outerSectorRowPitch();
  }
  int indexForRowPad(int row, int pad)
  {
    if (pad > numberOfPadsAtRow(row)) return -1;

    int index = 0;

    if (row > 0 && row <= numberOfInnerRows() )             for (int i = 1; i < row; i++) index += numberOfPadsAtRow(i);
    else if (row > numberOfInnerRows() && row <= numberOfRows()) for (int i = numberOfInnerRows() + 1; i < row; i++)  index += numberOfPadsAtRow(i);

    index += pad - 1;
    return index;
  }
};

struct St_tpcPadrowT0C : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcPadrowT0C, tpcPadrowT0> {
  float T0(int sector, int row) {return Struct(sector-1)->T0[row-1];}
};

struct St_TpcPhiDirectionC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcPhiDirectionC, tpcCorrection> {};

struct St_tpcPressureBC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_tpcPressureBC, tpcCorrection> {};

struct St_TpcrChargeC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcrChargeC, tpcCorrection> {};

struct St_tpcRDOMapC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOMapC, tpcRDOMap> {
  unsigned char 	nrows(int i = 0) 	const {return Struct(i)->nrows;}
  unsigned char 	index(int i = 0) 	const {return Struct(i)->idx;}
  unsigned char 	row(int i = 0) 	const {return Struct(i)->row;}
  unsigned char 	padMin(int i = 0) 	const {return Struct(i)->padMin;}
  unsigned char 	padMax(int i = 0) 	const {return Struct(i)->padMax;}
  unsigned char 	rdoI(int i = 0) 	const {return Struct(i)->rdo;}
  int         rdo(int padrow, int pad = 0) const;
};

struct St_tpcRDOMasksC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOMasksC, tpcRDOMasks>
{
  unsigned int 	runNumber(int i = 0) 	        {return Struct(i)->runNumber;}
  unsigned int 	sector(int i = 0) 	        {return Struct(i)->sector;}
  unsigned int 	mask(int i = 0) 	        {return Struct(i)->mask;}
  unsigned int        getSectorMask(unsigned int sector);
  static unsigned int rdoForPadrow(int row)   //Function returns the rdo board number for a given padrow index. Range of map used is 1-45.
  {
    unsigned int rdo = 0;

    if      (row > 0 && row <=  8) rdo = 1;
    else if (row > 8 && row <= 13) rdo = 2;
    else if (row > 13 && row <= 21) rdo = 3;
    else if (row > 21 && row <= 29) rdo = 4;
    else if (row > 29 && row <= 37) rdo = 5;
    else if (row > 37 && row <= 45) rdo = 6;

    return rdo;
  }
  static unsigned int rdoForPadrow(int sector, int row)   //Function returns the rdo board number for a given padrow index. Range of map used is 1-45.
  {
    if (St_tpcPadConfigC::instance()->iTpc(sector)) return 8;

    return rdoForPadrow(row);
  }
  bool        isOn(int sector, int rdo)
  {
    if (St_tpcPadConfigC::instance()->iTpc(sector)) return 1;

    if (sector < 1 || sector > 24 || rdo < 1 || rdo > 6)	return 0;

    unsigned int MASK = getSectorMask(sector);
    MASK = MASK >> (rdo - 1);
    MASK &= 0x00000001;
    return MASK;
  }
  bool       isRowOn(int sector, int row) {return isOn(sector, rdoForPadrow(sector, row));}
};

struct St_tpcRDOT0offsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcRDOT0offsetC, tpcRDOT0offset> {
  unsigned char* 	isShifted(int i = 0) 	const {return Struct(i)->isShifted;}
  bool        IsShfited(int sector) const {return isShifted()[sector-1];}
  float       T0(int sector, int padrow, int pad) const;
};

struct St_TpcRowQC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcRowQC, tpcCorrection> {};

struct St_tpcSCGLC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcSCGLC, tpcSCGL> {
  float* SC()                             {return Struct()->SC;}
  float* SCoffset()                       {return Struct()->SCoffset;}
  float* SCexponent()                     {return Struct()->SCexponent;}
  float* SCscaler()                       {return Struct()->SCscaler;}
  float* GL()                             {return Struct()->GL;}
  float* GLoffset()                       {return Struct()->GLoffset;}
  float  GLradius()                       {return Struct()->GLradius;}
  float  GLwidth()                        {return Struct()->GLwidth;}
  int    mode()                           {return Struct()->mode;}
  char*  comment()                        {return Struct()->comment;}
};

struct St_TpcSecRowCorC : tpcrs::IConfigStruct
{
  virtual TpcSecRowCor* Struct(int i = 0) const = 0;

  float* 	GainScale(int i = 0) 	        {return Struct(i)->GainScale;}
  float* 	GainRms(int i = 0) 	        {return Struct(i)->GainRms;}
};

struct St_TpcSecRowBC : tpcrs::ConfigStruct<St_TpcSecRowCorC, St_TpcSecRowBC, TpcSecRowCor> {};

struct St_TpcSecRowCC : tpcrs::ConfigStruct<St_TpcSecRowCorC, St_TpcSecRowCC, TpcSecRowCor> {};

struct St_tpcSectorT0offsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tpcSectorT0offsetC, tpcSectorT0offset> {
  float* 	t0(int i = 0) 	        {return Struct(i)->t0;}
  float       t0offset(int sector=1)        {return t0()[sector-1];}
};

struct St_TpcSpaceChargeC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcSpaceChargeC, tpcCorrection> {};

struct StTpcInnerSectorPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcInnerSectorPosition, Survey> {// Inner part of sector to Super Sector
};

struct StTpcOuterSectorPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcOuterSectorPosition, Survey> {// Outer part of sector to Super Sector
};

struct StTpcSuperSectorPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcSuperSectorPosition, Survey> {// Extra rotation for whole Super Sector to half Tpc
};

struct StTpcHalfPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcHalfPosition, Survey> {// Extra rotation for half of Tpc  to Tpc
  const TGeoHMatrix  &GetEastMatrix() {return  GetMatrix(TPC::Half::first);}
  const TGeoHMatrix  &GetWestMatrix() {return  GetMatrix(TPC::Half::second);}
  const TGeoHMatrix  &GetEastMatrixR() {return  GetMatrixR(TPC::Half::first);}
  const TGeoHMatrix  &GetWestMatrixR() {return  GetMatrixR(TPC::Half::second);}
  static void Normalize(TGeoHMatrix &R) {}
};

struct StTpcPosition : tpcrs::ConfigStruct<St_SurveyC, StTpcPosition, Survey> {// Global position of TPC in Magnet
  const TGeoHMatrix  &GetMatrix() {return  St_SurveyC::GetMatrix(0);}
};

struct St_TpcTanLC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcTanLC, tpcCorrection> {};

struct St_tpcTimeDependenceC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_tpcTimeDependenceC, tpcCorrection> {};

struct St_tpcWaterOutC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_tpcWaterOutC, tpcCorrection> {};

struct St_TpcZCorrectionBC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcZCorrectionBC, tpcCorrection> {};

struct St_TpcZDCC : tpcrs::ConfigStruct<St_tpcCorrectionC, St_TpcZDCC, tpcCorrection> {};

struct St_trgTimeOffsetC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_trgTimeOffsetC, trgTimeOffset>
{
  float 	offset(int i = 0)     	   {return Struct(i)->offset;}
  float 	laserOffset(int i = 0) 	   {return Struct(i)->laserOffset;}
  float 	laserOffsetW(int i = 0) 	   {return Struct(i)->laserOffsetW;}
  float       triggerTimeOffset(int i = 0)     {return 1e-6 * (mLaser ? laserOffset(i)  : offset(i));} // usec
  float       triggerTimeOffsetWest(int i = 0) {return 1e-6 * (mLaser ? laserOffsetW(i) :         0);} // usec
  void          SetLaser(bool k = kTRUE)         {mLaser = k;}
 private:
  bool        mLaser;
};

struct St_trigDetSumsC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_trigDetSumsC, trigDetSums>
{
  unsigned int 	runNumber(int i = 0) 	        {return Struct(i)->runNumber;}
  unsigned int 	timeOffset(int i = 0) 	{return Struct(i)->timeOffset;}
  double 	ctbWest(int i = 0) 	        {return Struct(i)->ctbWest;}
  double 	ctbEast(int i = 0) 	        {return Struct(i)->ctbEast;}
  double 	ctbTOFp(int i = 0) 	        {return Struct(i)->ctbTOFp;}
  double 	tofp(int i = 0) 	        {return Struct(i)->tofp;}
  double 	zdcWest(int i = 0) 	        {return Struct(i)->zdcWest;}
  double 	zdcEast(int i = 0) 	        {return Struct(i)->zdcEast;}
  double 	zdcX(int i = 0) 	        {return Struct(i)->zdcX;}
  double 	mult(int i = 0) 	        {return Struct(i)->mult;}
  double 	L0(int i = 0) 	        {return Struct(i)->L0;}
  double 	bbcX(int i = 0) 	        {return Struct(i)->bbcX;}
  double 	bbcXctbTOFp(int i = 0) 	{return Struct(i)->bbcXctbTOFp;}
  double 	bbcWest(int i = 0) 	        {return Struct(i)->bbcWest;}
  double 	bbcEast(int i = 0) 	        {return Struct(i)->bbcEast;}
  double 	bbcYellowBkg(int i = 0) 	{return Struct(i)->bbcYellowBkg;}
  double 	bbcBlueBkg(int i = 0) 	{return Struct(i)->bbcBlueBkg;}
  double 	pvpdWest(int i = 0) 	        {return Struct(i)->pvpdWest;}
  double 	pvpdEast(int i = 0) 	        {return Struct(i)->pvpdEast;}
  double 	zdcCoin(int i = 0)            {return Nc(zdcX(i), zdcEast(i), zdcWest(i));}
  double 	bbcCoin(int i = 0)            {return Nc(bbcX(i), bbcEast(i), bbcWest(i));}
  void		validityMargin(double margin = 0) {fMargin = margin;}
  double getCTBWest() {return ctbWest();}
  double getCTBEast() {return ctbEast();}
  double getCTBOrTOFp() {return ctbTOFp();}
  double getTOFp() {return tofp();}
  double getZDCWest() {return zdcWest();}
  double getZDCEast() {return zdcEast();}
  double getZDCX() {return zdcX();}
  double getZDCCoin() {return zdcCoin();}
  double getMult() {return mult();}
  double getL0() {return L0();}
  double getBBCX() {return bbcX();}
  double getBBCCoin() {return bbcCoin();}
  double getBBCXCTB() {return bbcXctbTOFp();}
  double getBBCWest() {return bbcWest();}
  double getBBCEast() {return bbcEast();}
  double getBBCYellowBkg() {return bbcYellowBkg();}
  double getBBCBlueBkg() {return bbcBlueBkg();}
  double getPVPDWest() {return pvpdWest();}
  double getPVPDEast() {return pvpdEast();}
  unsigned int   getRichHVStatus() {return St_richvoltagesC::instance()->status();}
  void     setValidityMargin(double margin = 0) {validityMargin(margin);}

  // The following code attempts to correct coincidence rates for accidentals and multiples
  // See STAR Note 528
  static double Nc(double New, double Ne, double Nw, int n_bunches = 111)
  {
    // 111 is a guess using the maximum seen filled bunches in RHIC so far
    // (not always the case, but we don't have access to this number)
    double Nbc = tpcrs::Cfg<starClockOnl>().frequency * ((double) n_bunches) / 120.;
    return -Nbc * TMath::Log(1. - ((New - (Ne * Nw / Nbc)) / (Nbc + New - Ne - Nw)));
  }
 private:
  double	fMargin;
};

float GainCorrection(int sector, int row);

#endif
