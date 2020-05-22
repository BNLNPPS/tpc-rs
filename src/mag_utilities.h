/***********************************************************************
 * Author: Jim Thomas   11/1/2000
 *
 * Description: Utilities for the Magnetic Field
 ***********************************************************************/
#ifndef TPCRS_MAG_UTILITIES_H_
#define TPCRS_MAG_UTILITIES_H_

#include <cmath>

#include "TArrayF.h"
#include "TArrayD.h"
#include "TMatrix.h"
#include "mag_field.h"

class TFile;
class TNtuple;

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
};

enum   CorrectSelect {
  kIterateUndo       = 0x1       // Bit 1
};

enum   EBMapSizes {
  BMap_nZ    =          57,           // Number of Z points in table. Measured STAR B field Maps from Steve T.
  BMap_nR    =          28,           // Number of R points in table.
  BMap_nPhi  =          37,           // Number of Phi points in table.
  EMap_nZ    =         224,           // Number of Z points in table. Standard STAR distortion tables for interpolating.
  EMap_nR    =          82,           // Number of R points in table
  EMap_nPhi  =          13            // Number of Phi points in table ( add one for 360 == 0 )
};

// DO NOT change the numbering of these constants. StBFChain depends
// on these values to build an option flag. The option flag used in
// the chain is 2x larger than shown here in order to allow the first
// bit to be used as an on/off flag.  It is shifted away before entering
// StMagUtilities.  So, this can be summarized by saying:
// Bit counting starts at 0 for the chain option flag (...,3,2,1,0)

class TDataSet ;
class St_tpcHVPlanesC;
class St_tpcCalibResolutionsC;
class St_tpcHighVoltagesC;
class St_tpcOmegaTauC;
class St_tpcGridLeakC;
class St_spaceChargeCorC;
class St_tpcChargeEventC;
class TRandom;


class StMagUtilities
{


 private:
  const tpcrs::Configurator& cfg_;
  const CoordTransform& transform_;
  StarMagField mag_field_;
  St_spaceChargeCorC*        fSpaceCharge   ;
  St_spaceChargeCorC*        fSpaceChargeR2 ;
  St_tpcHighVoltagesC*       fTpcVolts      ;
  St_tpcOmegaTauC*           fOmegaTau      ;
  St_tpcGridLeakC*           fGridLeak      ;
  St_tpcHVPlanesC*           fHVPlanes      ;
  St_tpcCalibResolutionsC*   fCalibResolutions ;
  St_tpcChargeEventC*        fAbortGapCharge;

  void    GetDistoSmearing ( int mode) ;
  void    GetMagFactor ()     ;
  void    GetTPCParams ()     ;
  void    GetTPCVoltages ( int mode ) ;
  void    GetSpaceCharge ()   ;
  void    GetSpaceChargeR2 () ;
  void    GetShortedRing ()   ;
  void    GetOmegaTau ()      ;
  void    GetGridLeak ( int mode ) ;
  void    GetHVPlanes()       ;
  void    GetE()              ;
  void    GetAbortGapCharge() ;

  void    CommonStart ( int mode ) ;
  int   IsPowerOfTwo (int i) ;
  void    SectorNumber ( int &Sector, const float x[] ) ;
  void    SectorNumber ( int &Sector, float phi, const float z ) ;
  void    GetGLWallData( const int select, float DataInTheGap[] ) ;
  int   SectorSide   ( int &Sector, const float x[] ) ;  // -1 for east, +1 for west
  int   SectorSide   ( int &Sector, const float z   ) ;
  float LimitZ (int &Sector, const float x[] ) ;
  float Interpolate2DTable  ( const int ORDER, const float x, const float y, const int nx, const int ny,
                                        const float XV[], const float YV[], const TMatrix &Array ) ;
  float Interpolate3DTable ( const int ORDER, const float x,    const float y,    const float z,
                                       const int  nx,    const int  ny,    const int  nz,
                                       const float XV[], const float YV[], const float ZV[],
                                       TMatrix** ArrayofArrays ) ;
  void    Interpolate2DEdistortion ( const int ORDER, const float r, const float z,
      const float Er[EMap_nZ][EMap_nR], float &Er_value ) ;
  void    Interpolate3DEdistortion ( const int ORDER, const float r, const float phi, const float z,
      const float Er[EMap_nZ][EMap_nPhi][EMap_nR], const float Ephi[EMap_nZ][EMap_nPhi][EMap_nR],
      float &Er_value, float &Ephi_value ) ;
  void    PoissonRelaxation  ( TMatrix &ArrayV, TMatrix &Charge, TMatrix &EroverEz,
                                       const int ITERATIONS ) ;

  void    Poisson3DRelaxation( TMatrix** ArrayofArrayV, TMatrix** ArrayofCharge, TMatrix** ArrayofEroverEz,
                                       TMatrix** ArrayofEPhioverEz,
                                       const int PHISLICES, const float DeltaPhi,
                                       const int ITERATIONS, const int SYMMETRY) ;

  int    mDistortionMode;             // Distortion mode - determines which corrections are run
  unsigned int   mCorrectionsMode;            // Corrections mode - determines how corrections are run
  bool   DoOnce ;                     // First pass: initializations

  float  StarDriftV ;                 // Drift Velocity (cm/microSec) Magnitude
  float  TPC_Z0 ;                     // Z location of STAR TPC Ground Wire Plane (cm) Magnitude
  float  XTWIST ;                     // X Displacement of West end of TPC wrt magnet (mRad)
  float  YTWIST ;                     // Y Displacement of West end of TPC wrt magnet (mRad)
  double CathodeV ;                   // Cathode Potential (volts)
  double GG ;                         // Gating Grid voltage (volts)
  float  GGideal ;                    // Ideal set GG voltage, not effective voltage
  float  Inner_GLW_Voltage[24] ;      // Voltage on the inside of the Grid Leak Wall facing the GG (~GG effective voltage)
  float  Outer_GLW_Voltage[24] ;      // Voltage on the outside surface of the Grid Leak Wall facing the outer sector
  float  Rtot ;                       // Total resistance of the (normal) resistor chain
  float  Rfrac ;                      // Fraction of full resistor chain inside TPC drift volume (~1.0)
  float  RPitch ;                     // Field Cage Ring to Ring pitch (cm)
  float  GGeffectiveness ;            // Effectiveness of GG voltage to be the average at its plane
  float  deltaGGeffectiveness ;       // Effectiveness of GG voltage changes to be expressed in average
  float  EASTCLOCKERROR ;             // Phi rotation of East end of TPC in milli-radians
  float  WESTCLOCKERROR ;             // Phi rotation of West end of TPC in milli-radians
  float  IFCRadius ;                  // Radius of the Inner Field Cage
  float  OFCRadius ;                  // Radius of the Outer Field Cage
  float  INNERGGFirst ;               // Radius of the first Inner Gating Grid Wire
  float  INNERGGLast ;                // Radius of the last Inner Gating Grid Wire
  float  OUTERGGFirst ;               // Radius of the first Outer Gating Grid Wire
  float  OUTERGGLast ;                // Radius of the last Outer Gating Grid Wire
  float  GAPRADIUS ;                  // Radius of the gap between the inner and outer grids (cm)
  float  WIREGAP ;                    // Width of the gap between the inner and outer grids (cm)
  double TPCROWR[24][128] ;           // Radii of TPC rows along the sector centerlines
  int    INNER[24];                   // Number of TPC rows in the inner sectors
  int    TPCROWS[24];                 // Total number of TPC rows per sector (Inner + Outer)
  float  StarMagE ;                   // STAR Electric Field (V/cm) Magnitude
  float  IFCShift ;                   // Shift of the IFC towards the West Endcap (cm)
  float  TensorV1 ;                   // Omega Tau tensor parameter - in the ExB direction
  float  TensorV2 ;                   // Omega Tau tensor parameter - in the direction perpendicular to ExB and Z axis
  float  Const_0, Const_1, Const_2  ; // OmegaTau parameters
  float  SpaceChargeEWRatio         ; // Ratio of East/West Space charge ... for example, d-Au should be ratio 6/5, Au-Au ratio 1/1
  double SpaceCharge, SpaceChargeR2 ; // Space Charge parameters (uniform or 1/R**2 in the TPC - arbitrary units)
  double InnerGridLeakStrength      ; // Relative strength of the Inner grid leak
  double InnerGridLeakRadius        ; // Location (in local Y coordinates) of the Inner grid leak
  double InnerGridLeakWidth         ; // Half-width of the Inner grid leak.  Must be larger than life for numerical reasons.
  double MiddlGridLeakStrength      ; // Relative strength of the Middle grid leak
  double MiddlGridLeakRadius        ; // Location (in local Y coordinates) of the Middle grid leak
  double MiddlGridLeakWidth         ; // Half-width of the Middle grid leak.  Must be larger than life for numerical reasons.
  double OuterGridLeakStrength      ; // Relative strength of the Outer grid leak
  double OuterGridLeakRadius        ; // Location (in local Y coordinates) of the Outer grid leak
  double OuterGridLeakWidth         ; // Half-width of the Outer grid leak.  Must be larger than life for numerical reasons.
  float  GLWeights[96]              ; // GridLeak weights per sector.  24 sectors x 3 locations
  int    ShortTableRows             ; // Number of rows in the Shorted Ring Table
  int    Side[10]                   ; // Location of Short   E=0 /   W=1
  int    Cage[10]                   ; // Location of Short IFC=0 / OFC=1
  float  Ring[10]                   ; // Location of Short counting out from the CM.  CM==0
  float  MissingResistance[10]      ; // Amount of Missing Resistance due to this short (MOhm)
  float  Resistor[10]               ; // Amount of compensating resistance added for this short
  float  deltaVGGEast               ; // Voltage error on the East Gated Grid
  float  deltaVGGWest               ; // Voltage error on the West Gated Grid
  bool   iterateDistortion          ; // Flag on whether to iterate in determining distortions
  int    iterationFailCounter       ; // Count of number of iteration fails
  bool   doingDistortion            ; // Flag on whether doing or undoing distortions
  bool   usingCartesian             ; // Using Cartesian or cylindrical coordinates
  TRandom* mRandom                    ; // Random number generator (used in distortion smearing)
  float  SmearCoefSC                ; // Distortion smearing coefficient for SpaceCharge
  float  SmearCoefGL                ; // Distortion smearing coefficient for GridLeak
  TArrayF* AbortGapCharges            ; // Charges deposited into the TPC due to Abort Gap Cleaning events
  TArrayD* AbortGapTimes              ; // Times since charges deposited into the TPC due to Abort Gap Cleaning events
  float  AbortGapChargeCoef         ; // Scale factor for charge deposited due to Abort Gap Cleaning events
  float  IonDriftVel                ; // Drift velocity of ions in the TPC gas



  float  shiftEr[EMap_nZ][EMap_nR] ;
  float  spaceEr[EMap_nZ][EMap_nR] ;
  float  spaceR2Er[EMap_nZ][EMap_nR] ;
  float  shortEr[EMap_nZ][EMap_nR] ;
  float  GGVoltErrorEr[EMap_nZ][EMap_nR] ;

  static   float ePhiList[EMap_nPhi] ;   // Note: These are initialized near CommonStart() in the .cxx file
  static   float eRList[EMap_nR]     ;
  static   float eZList[EMap_nZ]     ;
  static   TNtuple* fgDoDistortion;
  static   TNtuple* fgUnDoDistortion;
 public:

  StMagUtilities(const tpcrs::Configurator& cfg, const CoordTransform& trans, int mode = 0);
  StMagUtilities(const tpcrs::Configurator& cfg, const CoordTransform& trans, const StarMagField::EBField map, const float factor, int mode);

  void    B3DFieldTpc ( const float xTpc[], float BTpc[], int Sector = -1 );
  void    BFieldTpc ( const float xTpc[], float BTpc[], int Sector = -1 );

  void    DoDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoBDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    Undo2DBDistortion ( const float x[], float Xprime[], int Sector = -1 ); // {UndoBDistortion(x,Xprime,Sector);}
  void    FastUndoBDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    FastUndo2DBDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoPad13Distortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoPad40Distortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoTwistDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoClockDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoMembraneDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoEndcapDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoSpaceChargeDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoSpaceChargeR0Distortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoSpaceChargeR2Distortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoGridLeakDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    Undo2DGridLeakDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    Undo3DGridLeakDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoFullGridLeakDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoIFCShiftDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoShortedRingDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoGGVoltErrorDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoSectorAlignDistortion ( const float x[], float Xprime[], int Sector = -1 ) ;
  void    UndoAbortGapDistortion ( const float x[], float Xprime[], int Sector = -1, float TimeSinceDeposition = -1.0 ) ;

  void    FixSpaceChargeDistortion ( const int Charge, const float x[3], const float p[3],
      const Prime PrimaryOrGlobal,
      float x_new[3], float p_new[3],
      const unsigned int RowMask1 = 0xFFFFFF00,
      const unsigned int RowMask2 = 0x1FFFFF,
      const float VertexError = 0.0200 ) ;

  void    ApplySpaceChargeDistortion ( const double sc, const int Charge,
      const float x[3], const float p[3],
      const Prime PrimaryOrGlobal, int &new_Charge,
      float x_new[3], float p_new[3],
      const unsigned int RowMask1 = 0xFFFFFF00,
      const unsigned int RowMask2 = 0x1FFFFF,
      const float VertexError = 0.0200 ) ;

  int   PredictSpaceChargeDistortion ( int   sec,
      int   Charge,
      float Pt,
      float VertexZ,
      float PseudoRapidity,
      float DCA,
      const unsigned int RowMask1,
      const unsigned int RowMask2,
      float &pSpace ) ;

  int   PredictSpaceChargeDistortion ( int   sec,
      int   Charge,
      float Pt,
      float VertexZ,
      float PseudoRapidity,
      float Phi,
      float DCA,
      const unsigned long long RowMask1,
      const unsigned long long RowMask2,
      float RowMaskErrorR[64],
      float RowMaskErrorRPhi[64],
      float &pSpace ) ;

  int   PredictSpaceChargeDistortion ( int   NHits,
      int   Charge,
      float Pt,
      float VertexZ,
      float PseudoRapidity,
      float Phi,
      float DCA,
      double R[128],
      double ErrorR[128],
      double ErrorRPhi[128],
      float &pSpace ) ;

  void     ManualSpaceChargeR2(double SpcChg, float EWRatio = 1.0 )
  {
    SpaceChargeR2 = SpcChg ; fSpaceChargeR2 = 0 ;
    SpaceChargeEWRatio = EWRatio ;
  }
  void     AutoSpaceCharge()   {GetSpaceCharge()  ; } // use DB
  void     AutoSpaceChargeR2() {GetSpaceChargeR2(); } // use DB
  double CurrentSpaceCharge()   {return SpaceCharge  ;}
  double CurrentSpaceChargeR2() {return SpaceChargeR2;}
  float  CurrentSpaceChargeEWRatio() { return SpaceChargeEWRatio ; }
  bool   UpdateTPCHighVoltages();
  bool   UpdateShortedRing();
  void     UseIterativeUndoDistortion(bool flag = true) { iterateDistortion = flag; }
  int    IterationFailCount(); // must be called once before first actual use
  float  GetConst_0() { return Const_0; }
  float  GetConst_1() { return Const_1; }
  float  GetConst_2() { return Const_2; }
  static  void    SetDoDistortionT  (TFile* f = 0);
  static  void    SetUnDoDistortionT(TFile* f = 0);

  void     Cart2Polar(const float* x, float &r, float &phi)
  {
    r      =  std::sqrt( x[0] * x[0] + x[1] * x[1] ) ;
    phi    =  std::atan2(x[1], x[0]) ;
  }
  void     Cart2Polar(const float* x, double &r, double &phi)
  {
    r      =  std::sqrt( x[0] * x[0] + x[1] * x[1] ) ;
    phi    =  std::atan2(x[1], x[0]) ;
  }
  void     Polar2Cart(const double r, const double phi, float* Xprime)
  {
    Xprime[0] = r * std::cos(phi) ;
    Xprime[1] = r * std::sin(phi) ;
  }
};

#endif
