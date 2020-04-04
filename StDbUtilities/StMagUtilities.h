/***********************************************************************
 * Author: Jim Thomas   11/1/2000
 *
 * Description: Utilities for the Magnetic Field
 ***********************************************************************/
#ifndef StMagUtilities_H
#define StMagUtilities_H
#define  __NEW_MagUtilities__
#include <cmath>

#include "TArrayF.h"
#include "TArrayD.h"
#include "TMatrix.h"      // TMatrix keeps changing ... keep it here until proven otherwise.
#include "StarMagField/StarMagField.h"
class TFile;
class TNtuple;
#include "TMath.h"

enum   EBField  { kUndefined = 0, kConstant = 1, kMapped = 2, kChain = 3 } ;
enum   Prime    { IsPrimary = 0, IsGlobal = 1 } ;

// Bit counting starts at 1 for the mode switch (...,3,2,1)

enum   DistortSelect {
  kBMap              = 0x08,     // Bit 4
  kPadrow13          = 0x10,     // Bit 5
  kTwist             = 0x20,     // Bit 6
  kClock             = 0x40,     // Bit 7
  kMembrane          = 0x80,     // Bit 8
  kEndcap            = 0x100,    // Bit 9
  kIFCShift          = 0x200,    // Bit 10
  kSpaceCharge       = 0x400,    // Bit 11
  kSpaceChargeR2     = 0x800,    // Bit 12
  kShortedRing       = 0x1000,   // Bit 13
  kFast2DBMap        = 0x2000,   // Bit 14
  kGridLeak          = 0x4000,   // Bit 15
  k3DGridLeak        = 0x8000,   // Bit 16
  kGGVoltError       = 0x10000,  // Bit 17
  kSectorAlign       = 0x20000,  // Bit 18
  kDisableTwistClock = 0x40000,  // Bit 19
  kFullGridLeak      = 0x80000,  // Bit 20
  kDistoSmearing     = 0x100000, // Bit 21
  kPadrow40          = 0x200000, // Bit 22
  kAbortGap          = 0x400000  // Bit 23
} ;
enum   CorrectSelect {
  kIterateUndo       = 0x1       // Bit 1
} ;
enum   EBMapSizes {
  BMap_nZ    =          57,           // Number of Z points in table. Measured STAR B field Maps from Steve T.
  BMap_nR    =          28,           // Number of R points in table.
  BMap_nPhi  =          37,           // Number of Phi points in table.
  EMap_nZ    =         224,           // Number of Z points in table. Standard STAR distortion tables for interpolating.
  EMap_nR    =          82,           // Number of R points in table
  EMap_nPhi  =          13            // Number of Phi points in table ( add one for 360 == 0 )
} ;

// DO NOT change the numbering of these constants. StBFChain depends
// on these values to build an option flag. The option flag used in
// the chain is 2x larger than shown here in order to allow the first
// bit to be used as an on/off flag.  It is shifted away before entering
// StMagUtilities.  So, this can be summarized by saying:
// Bit counting starts at 0 for the chain option flag (...,3,2,1,0)

class StTpcDb ;
class TDataSet ;
class St_tpcHVPlanesC;
class St_tpcCalibResolutionsC;
class St_tpcHighVoltagesC;
class St_tpcOmegaTauC;
class St_tpcGridLeakC;
class St_spaceChargeCorC;
class St_tpcChargeEventC;
class TRandom;

//class TMatrix ;

class StMagUtilities
{


 private:
  static StMagUtilities* fgInstance;
  St_spaceChargeCorC*        fSpaceCharge   ;
  St_spaceChargeCorC*        fSpaceChargeR2 ;
  St_tpcHighVoltagesC*       fTpcVolts      ;
  St_tpcOmegaTauC*           fOmegaTau      ;
  St_tpcGridLeakC*           fGridLeak      ;
  St_tpcHVPlanesC*           fHVPlanes      ;
  St_tpcCalibResolutionsC*   fCalibResolutions ;
  St_tpcChargeEventC*        fAbortGapCharge;

  virtual void    GetDistoSmearing ( Int_t mode) ;
  virtual void    GetMagFactor ()     ;
  virtual void    GetTPCParams ()     ;
  virtual void    GetTPCVoltages ( Int_t mode ) ;
  virtual void    GetSpaceCharge ()   ;
  virtual void    GetSpaceChargeR2 () ;
  virtual void    GetShortedRing ()   ;
  virtual void    GetOmegaTau ()      ;
  virtual void    GetGridLeak ( Int_t mode ) ;
  virtual void    GetHVPlanes()       ;
  virtual void    GetE()              ;
  virtual void    GetAbortGapCharge() ;

  virtual void    CommonStart ( Int_t mode ) ;
  virtual void    Search ( const Int_t N, const Float_t Xarray[], const Float_t x, Int_t &low )
  {StarMagField::Instance()->Search(N, Xarray, x, low);}
  virtual Int_t   IsPowerOfTwo (Int_t i) ;
  virtual void    SectorNumber ( Int_t &Sector, const Float_t x[] ) ;
  virtual void    SectorNumber ( Int_t &Sector, Float_t phi, const Float_t z ) ;
  virtual void    GetGLWallData( const Int_t select, Float_t DataInTheGap[] ) ;
  virtual Int_t   SectorSide   ( Int_t &Sector, const Float_t x[] ) ;  // -1 for east, +1 for west
  virtual Int_t   SectorSide   ( Int_t &Sector, const Float_t z   ) ;
  virtual Float_t LimitZ (Int_t &Sector, const Float_t x[] ) ;
  virtual Float_t Interpolate ( const Float_t Xarray[], const Float_t Yarray[],
                                const Int_t ORDER, const Float_t x )
  {return StarMagField::Instance()->Interpolate(Xarray, Yarray, ORDER, x);}
  virtual Float_t Interpolate2DTable  ( const Int_t ORDER, const Float_t x, const Float_t y, const Int_t nx, const Int_t ny,
                                        const Float_t XV[], const Float_t YV[], const TMatrix &Array ) ;
  virtual Float_t Interpolate3DTable ( const Int_t ORDER, const Float_t x,    const Float_t y,    const Float_t z,
                                       const Int_t  nx,    const Int_t  ny,    const Int_t  nz,
                                       const Float_t XV[], const Float_t YV[], const Float_t ZV[],
                                       TMatrix** ArrayofArrays ) ;
  virtual void    Interpolate2DBfield ( const Float_t r, const Float_t z,
                                        Float_t &Br_value, Float_t &Bz_value )
  {StarMagField::Instance()->Interpolate2DBfield ( r, z, Br_value, Bz_value );}
  virtual void    Interpolate3DBfield ( const Float_t r, const Float_t z, const Float_t phi,
                                        Float_t &Br_value, Float_t &Bz_value, Float_t &Bphi_value )
  {StarMagField::Instance()->Interpolate3DBfield ( r, z, phi, Br_value, Bz_value, Bphi_value );}
  virtual void    Interpolate2DEdistortion ( const Int_t ORDER, const Float_t r, const Float_t z,
      const Float_t Er[EMap_nZ][EMap_nR], Float_t &Er_value ) ;
  virtual void    Interpolate3DEdistortion ( const Int_t ORDER, const Float_t r, const Float_t phi, const Float_t z,
      const Float_t Er[EMap_nZ][EMap_nPhi][EMap_nR], const Float_t Ephi[EMap_nZ][EMap_nPhi][EMap_nR],
      Float_t &Er_value, Float_t &Ephi_value ) ;
  virtual void    PoissonRelaxation  ( TMatrix &ArrayV, TMatrix &Charge, TMatrix &EroverEz,
                                       const Int_t ITERATIONS ) ;

  virtual void    Poisson3DRelaxation( TMatrix** ArrayofArrayV, TMatrix** ArrayofCharge, TMatrix** ArrayofEroverEz,
                                       TMatrix** ArrayofEPhioverEz,
                                       const Int_t PHISLICES, const Float_t DeltaPhi,
                                       const Int_t ITERATIONS, const Int_t SYMMETRY) ;

  Int_t    mDistortionMode;             // Distortion mode - determines which corrections are run
  UInt_t   mCorrectionsMode;            // Corrections mode - determines how corrections are run
  Bool_t   DoOnce ;                     // First pass: initializations

  Float_t  StarDriftV ;                 // Drift Velocity (cm/microSec) Magnitude
  Float_t  TPC_Z0 ;                     // Z location of STAR TPC Ground Wire Plane (cm) Magnitude
  Float_t  XTWIST ;                     // X Displacement of West end of TPC wrt magnet (mRad)
  Float_t  YTWIST ;                     // Y Displacement of West end of TPC wrt magnet (mRad)
  Double_t CathodeV ;                   // Cathode Potential (volts)
  Double_t GG ;                         // Gating Grid voltage (volts)
  Float_t  GGideal ;                    // Ideal set GG voltage, not effective voltage
  Float_t  Inner_GLW_Voltage[24] ;      // Voltage on the inside of the Grid Leak Wall facing the GG (~GG effective voltage)
  Float_t  Outer_GLW_Voltage[24] ;      // Voltage on the outside surface of the Grid Leak Wall facing the outer sector
  Float_t  Rtot ;                       // Total resistance of the (normal) resistor chain
  Float_t  Rfrac ;                      // Fraction of full resistor chain inside TPC drift volume (~1.0)
  Float_t  RPitch ;                     // Field Cage Ring to Ring pitch (cm)
  Float_t  GGeffectiveness ;            // Effectiveness of GG voltage to be the average at its plane
  Float_t  deltaGGeffectiveness ;       // Effectiveness of GG voltage changes to be expressed in average
  Float_t  EASTCLOCKERROR ;             // Phi rotation of East end of TPC in milli-radians
  Float_t  WESTCLOCKERROR ;             // Phi rotation of West end of TPC in milli-radians
  Float_t  IFCRadius ;                  // Radius of the Inner Field Cage
  Float_t  OFCRadius ;                  // Radius of the Outer Field Cage
  Float_t  INNERGGFirst ;               // Radius of the first Inner Gating Grid Wire
  Float_t  INNERGGLast ;                // Radius of the last Inner Gating Grid Wire
  Float_t  OUTERGGFirst ;               // Radius of the first Outer Gating Grid Wire
  Float_t  OUTERGGLast ;                // Radius of the last Outer Gating Grid Wire
  Float_t  GAPRADIUS ;                  // Radius of the gap between the inner and outer grids (cm)
  Float_t  WIREGAP ;                    // Width of the gap between the inner and outer grids (cm)
  Double_t TPCROWR[24][128] ;           // Radii of TPC rows along the sector centerlines
  Int_t    INNER[24];                   // Number of TPC rows in the inner sectors
  Int_t    TPCROWS[24];                 // Total number of TPC rows per sector (Inner + Outer)
  Float_t  StarMagE ;                   // STAR Electric Field (V/cm) Magnitude
  Float_t  IFCShift ;                   // Shift of the IFC towards the West Endcap (cm)
  Float_t  TensorV1 ;                   // Omega Tau tensor parameter - in the ExB direction
  Float_t  TensorV2 ;                   // Omega Tau tensor parameter - in the direction perpendicular to ExB and Z axis
  Float_t  Const_0, Const_1, Const_2  ; // OmegaTau parameters
  Float_t  SpaceChargeEWRatio         ; // Ratio of East/West Space charge ... for example, d-Au should be ratio 6/5, Au-Au ratio 1/1
  Double_t SpaceCharge, SpaceChargeR2 ; // Space Charge parameters (uniform or 1/R**2 in the TPC - arbitrary units)
  Double_t InnerGridLeakStrength      ; // Relative strength of the Inner grid leak
  Double_t InnerGridLeakRadius        ; // Location (in local Y coordinates) of the Inner grid leak
  Double_t InnerGridLeakWidth         ; // Half-width of the Inner grid leak.  Must be larger than life for numerical reasons.
  Double_t MiddlGridLeakStrength      ; // Relative strength of the Middle grid leak
  Double_t MiddlGridLeakRadius        ; // Location (in local Y coordinates) of the Middle grid leak
  Double_t MiddlGridLeakWidth         ; // Half-width of the Middle grid leak.  Must be larger than life for numerical reasons.
  Double_t OuterGridLeakStrength      ; // Relative strength of the Outer grid leak
  Double_t OuterGridLeakRadius        ; // Location (in local Y coordinates) of the Outer grid leak
  Double_t OuterGridLeakWidth         ; // Half-width of the Outer grid leak.  Must be larger than life for numerical reasons.
  Float_t  GLWeights[96]              ; // GridLeak weights per sector.  24 sectors x 3 locations
  Int_t    ShortTableRows             ; // Number of rows in the Shorted Ring Table
  Int_t    Side[10]                   ; // Location of Short   E=0 /   W=1
  Int_t    Cage[10]                   ; // Location of Short IFC=0 / OFC=1
  Float_t  Ring[10]                   ; // Location of Short counting out from the CM.  CM==0
  Float_t  MissingResistance[10]      ; // Amount of Missing Resistance due to this short (MOhm)
  Float_t  Resistor[10]               ; // Amount of compensating resistance added for this short
  Float_t  deltaVGGEast               ; // Voltage error on the East Gated Grid
  Float_t  deltaVGGWest               ; // Voltage error on the West Gated Grid
  Bool_t   useManualSCForPredict      ; // Flag on using fixed SC value or manually set one for Predict()
  Bool_t   iterateDistortion          ; // Flag on whether to iterate in determining distortions
  Int_t    iterationFailCounter       ; // Count of number of iteration fails
  Bool_t   doingDistortion            ; // Flag on whether doing or undoing distortions
  Bool_t   usingCartesian             ; // Using Cartesian or cylindrical coordinates
  TRandom* mRandom                    ; // Random number generator (used in distortion smearing)
  Float_t  SmearCoefSC                ; // Distortion smearing coefficient for SpaceCharge
  Float_t  SmearCoefGL                ; // Distortion smearing coefficient for GridLeak
  TArrayF* AbortGapCharges            ; // Charges deposited into the TPC due to Abort Gap Cleaning events
  TArrayD* AbortGapTimes              ; // Times since charges deposited into the TPC due to Abort Gap Cleaning events
  Float_t  AbortGapChargeCoef         ; // Scale factor for charge deposited due to Abort Gap Cleaning events
  Float_t  IonDriftVel                ; // Drift velocity of ions in the TPC gas



  Float_t  shiftEr[EMap_nZ][EMap_nR] ;
  Float_t  spaceEr[EMap_nZ][EMap_nR] ;
  Float_t  spaceR2Er[EMap_nZ][EMap_nR] ;
  Float_t  shortEr[EMap_nZ][EMap_nR] ;
  Float_t  GGVoltErrorEr[EMap_nZ][EMap_nR] ;

  static   Float_t ePhiList[EMap_nPhi] ;   // Note: These are initialized near CommonStart() in the .cxx file
  static   Float_t eRList[EMap_nR]     ;
  static   Float_t eZList[EMap_nZ]     ;
  static   TNtuple* fgDoDistortion;
  static   TNtuple* fgUnDoDistortion;
 public:

  StMagUtilities ( StTpcDb* dbin = 0, Int_t mode = 0 ) ;
  StMagUtilities ( const StarMagField::EBField map, const Float_t factor, Int_t mode );
  virtual ~StMagUtilities () {fgInstance = 0;}
  static StMagUtilities* Instance();

  virtual void    BField ( const Float_t x[], Float_t B[] )
  {StarMagField::Instance()->BField(x, B);}
  virtual void    BrBzField( const Float_t r, const Float_t z, Float_t &Br_value, Float_t &Bz_value )
  {StarMagField::Instance()->BrBzField(r, z, Br_value, Bz_value );}
  virtual void    B3DField ( const Float_t x[], Float_t B[] )
  {StarMagField::Instance()->B3DField(x, B);}
  virtual void    B3DFieldTpc ( const Float_t xTpc[], Float_t BTpc[], Int_t Sector = -1 );
  virtual void    BFieldTpc ( const Float_t xTpc[], Float_t BTpc[], Int_t Sector = -1 );
  virtual void    BrBz3DField ( const Float_t r, const Float_t z, const Float_t phi,
                                Float_t &Br_value, Float_t &Bz_value, Float_t &Bphi_value )
  {StarMagField::Instance()->BrBz3DField(r, z, phi, Br_value, Bz_value, Bphi_value);}

  virtual bool    UsingDistortion( const DistortSelect distortion ) { return ((mDistortionMode & distortion) ? true : false); }

  virtual void    DoDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoBDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    Undo2DBDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ); // {UndoBDistortion(x,Xprime,Sector);}
  virtual void    FastUndoBDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    FastUndo2DBDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoPad13Distortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoPad40Distortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoTwistDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoClockDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoMembraneDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoEndcapDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoSpaceChargeDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoSpaceChargeR0Distortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoSpaceChargeR2Distortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoGridLeakDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    Undo2DGridLeakDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    Undo3DGridLeakDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoFullGridLeakDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoIFCShiftDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoShortedRingDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoGGVoltErrorDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoSectorAlignDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1 ) ;
  virtual void    UndoAbortGapDistortion ( const Float_t x[], Float_t Xprime[], Int_t Sector = -1, Float_t TimeSinceDeposition = -1.0 ) ;

  virtual void    FixSpaceChargeDistortion ( const Int_t Charge, const Float_t x[3], const Float_t p[3],
      const Prime PrimaryOrGlobal,
      Float_t x_new[3], Float_t p_new[3],
      const unsigned int RowMask1 = 0xFFFFFF00,
      const unsigned int RowMask2 = 0x1FFFFF,
      const Float_t VertexError = 0.0200 ) ;

  virtual void    ApplySpaceChargeDistortion ( const Double_t sc, const Int_t Charge,
      const Float_t x[3], const Float_t p[3],
      const Prime PrimaryOrGlobal, Int_t &new_Charge,
      Float_t x_new[3], Float_t p_new[3],
      const unsigned int RowMask1 = 0xFFFFFF00,
      const unsigned int RowMask2 = 0x1FFFFF,
      const Float_t VertexError = 0.0200 ) ;

  virtual Int_t   PredictSpaceChargeDistortion ( Int_t   sec,
      Int_t   Charge,
      Float_t Pt,
      Float_t VertexZ,
      Float_t PseudoRapidity,
      Float_t DCA,
      const unsigned int RowMask1,
      const unsigned int RowMask2,
      Float_t &pSpace ) ;

  virtual Int_t   PredictSpaceChargeDistortion ( Int_t   sec,
      Int_t   Charge,
      Float_t Pt,
      Float_t VertexZ,
      Float_t PseudoRapidity,
      Float_t Phi,
      Float_t DCA,
      const unsigned long long RowMask1,
      const unsigned long long RowMask2,
      Float_t RowMaskErrorR[64],
      Float_t RowMaskErrorRPhi[64],
      Float_t &pSpace ) ;

  virtual Int_t   PredictSpaceChargeDistortion ( Int_t   NHits,
      Int_t   Charge,
      Float_t Pt,
      Float_t VertexZ,
      Float_t PseudoRapidity,
      Float_t Phi,
      Float_t DCA,
      Double_t R[128],
      Double_t ErrorR[128],
      Double_t ErrorRPhi[128],
      Float_t &pSpace ) ;

  virtual void     ManualShortedRing ( Int_t EastWest, Int_t InnerOuter,
                                       Float_t RingNumber, Float_t MissingRValue, Float_t ExtraRValue) ;

  virtual Int_t    GetSpaceChargeMode();
  virtual void     ManualSpaceCharge(Double_t SpcChg)   { SpaceCharge   = SpcChg ; fSpaceCharge   = 0 ; }
  virtual void     ManualSpaceChargeR2(Double_t SpcChg, Float_t EWRatio = 1.0 )
  {
    SpaceChargeR2 = SpcChg ; fSpaceChargeR2 = 0 ;
    SpaceChargeEWRatio = EWRatio ;
  }
  virtual void     ManualGridLeakStrength(Double_t inner, Double_t middle, Double_t outer);
  virtual void     ManualGridLeakRadius  (Double_t inner, Double_t middle, Double_t outer);
  virtual void     ManualGridLeakWidth   (Double_t inner, Double_t middle, Double_t outer);
  virtual void     AutoSpaceCharge()   {GetSpaceCharge()  ; } // use DB
  virtual void     AutoSpaceChargeR2() {GetSpaceChargeR2(); } // use DB
  virtual Double_t CurrentSpaceCharge()   {return SpaceCharge  ;}
  virtual Double_t CurrentSpaceChargeR2() {return SpaceChargeR2;}
  virtual Float_t  CurrentSpaceChargeEWRatio() { return SpaceChargeEWRatio ; }
  virtual Bool_t   UpdateTPCHighVoltages();
  virtual Bool_t   UpdateShortedRing();
  virtual void     UseManualSCForPredict(Bool_t flag = kTRUE) { useManualSCForPredict = flag; }
  virtual void     ManualGGVoltError(Double_t east, Double_t west);
  virtual void     UseIterativeUndoDistortion(Bool_t flag = kTRUE) { iterateDistortion = flag; }
  virtual Int_t    IterationFailCount(); // must be called once before first actual use
  Float_t  GetConst_0() { return Const_0; }
  Float_t  GetConst_1() { return Const_1; }
  Float_t  GetConst_2() { return Const_2; }
  static  void    SetDoDistortionT  (TFile* f = 0);
  static  void    SetUnDoDistortionT(TFile* f = 0);

  virtual void     Cart2Polar(const Float_t* x, Float_t &r, Float_t &phi)
  {
    r      =  TMath::Sqrt( x[0] * x[0] + x[1] * x[1] ) ;
    phi    =  TMath::ATan2(x[1], x[0]) ;
  }
  virtual void     Cart2Polar(const Float_t* x, Double_t &r, Double_t &phi)
  {
    r      =  TMath::Sqrt( x[0] * x[0] + x[1] * x[1] ) ;
    phi    =  TMath::ATan2(x[1], x[0]) ;
  }
  virtual void     Polar2Cart(const Double_t r, const Double_t phi, Float_t* Xprime)
  {
    Xprime[0] = r * TMath::Cos(phi) ;
    Xprime[1] = r * TMath::Sin(phi) ;
  }




};

#endif









