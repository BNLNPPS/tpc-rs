/***********************************************************************
 * Author: Jim Thomas   11/1/2000
 *
 * Description: the STAR Magnetic Field
  ***********************************************************************/

#ifndef StarMagField_H
#define StarMagField_H

#include <cmath>
#include <Rtypes.h>
#include "TH2.h"
#include "TGeoMatrix.h"

class StarMagField
{
 public:
  enum   EBField  { kUndefined = 0, kConstant = 1, kMapped = 2, kChain = 3 } ;
  enum   ESmFSizes {nZ = 57, nR = 28, nPhi = 37, nZSteel = 16, nRSteel = 115, nPhiSteel = 25};
  static  void    Search ( Int_t N, const Float_t Xarray[], Float_t x, Int_t &low ) ;
  virtual Float_t Interpolate ( const Float_t Xarray[], const Float_t Yarray[],
                                const Int_t ORDER, const Float_t x ) ;
  virtual void    Interpolate2DBfield ( const Float_t r, const Float_t z,
                                        Float_t &Br_value, Float_t &Bz_value ) ;
  virtual void    Interpolate2ExtDBfield ( const Float_t r, const Float_t z,
      Float_t &Br_value, Float_t &Bz_value ) ;
  virtual void    Interpolate3DBfield ( const Float_t r, const Float_t z, const Float_t phi,
                                        Float_t &Br_value, Float_t &Bz_value, Float_t &Bphi_value ) ;

 private:

  virtual void    ReadField ( ) ;
  static StarMagField* fgInstance;
  TGeoRotation fStarMagFieldRotation;
  TH2F* fBzdZCorrection; // correction due to endcap calomiter
  TH2F* fBrdZCorrection; // correction due to endcap calomiter
 public:
  //added by Lijuan

  virtual void    Interpolate3DBSteelfield ( const Float_t r, const Float_t z, const Float_t phi,
      Float_t &Br_value, Float_t &Bz_value, Float_t &Bphi_value );
  //end added by Lijuan



  EBField  fMap;       // (D) = kMapped; Global flag to indicate static arrays are full
  Float_t  fFactor;    // (D) = 1.0    ; Multiplicative factor (allows scaling and sign reversal)
  Float_t  fRescale;   // (D) = 1.0    ; Multiplicative factor (allows re-scaling wrt which map read)
  Float_t  fBDipole;   // (D) = -42.67 ; field value (kG)
  Float_t  fRmaxDip;   // (D) =  15.34 ; Inner field volume radius
  Float_t  fZminDip;   // (D) =  980.0 ; StArt of the DX mAgnet in Z
  Float_t  fZmaxDip;   // (D) = 1350.0 ; End of the DX mAgnet in Z
  Bool_t   fLock;      // (D) = kFALSE ; Set kTRUE if lock above values

  Float_t  Bz[nZ][nR], Br[nZ][nR] ;
  Float_t  Radius[nR], ZList[nZ] ;
  Float_t  Bz3D[nPhi][nZ][nR], Br3D[nPhi][nZ][nR], Bphi3D[nPhi][nZ][nR] ;
  Float_t  R3D[nR], Z3D[nZ], Phi3D[nPhi] ;
  Float_t  R3DSteel[nRSteel], Z3DSteel[nZSteel], Phi3DSteel[nPhiSteel] ;
  Float_t  Bz3DSteel[nPhiSteel][nZSteel][nRSteel];
  Float_t  Bx3DSteel[nPhiSteel][nZSteel][nRSteel], By3DSteel[nPhiSteel][nZSteel][nRSteel] ;
  static bool     mConstBz;

  //added by Lijuan
  Float_t  Br3DSteel[nPhiSteel][nZSteel][nRSteel], Bphi3DSteel[nPhiSteel][nZSteel][nRSteel] ;
  //end added by Lijuan

 public:

  StarMagField ( EBField map     = kMapped, Float_t Factor  =      1,
                 Bool_t  Lock    =  kFALSE, Float_t Rescale =      1,
                 Float_t Bdipole =  -42.67, Float_t Rmaxdip =  15.34,
                 Float_t Zmindip =   980.0, Float_t Zmaxdip = 1350.0) ;
  virtual ~StarMagField ()
  {
    fgInstance = 0;
    SafeDelete(fBzdZCorrection);
    SafeDelete(fBrdZCorrection);
  }
  static StarMagField* Instance();

  static void setConstBz( bool state ) { mConstBz = state; }

  virtual void    BField   ( const Float_t x[], Float_t B[] ) ;
  virtual void    BField   ( const Double_t x[], Double_t B[] ) ;
  virtual void    Field    ( const Float_t x[], Float_t B[] ) {BField(x, B);}
  virtual void    Field    ( const Double_t x[], Double_t B[] ) {BField(x, B);}
  virtual void    BrBzField( const Float_t r, const Float_t z, Float_t &Br_value, Float_t &Bz_value ) ;
  virtual void    B3DField ( const Float_t x[], Float_t B[] ) ;
  virtual void    B3DField ( const Double_t x[], Double_t B[] ) ;
  virtual void    BrBz3DField ( const Float_t r, const Float_t z, const Float_t phi,
                                Float_t &Br_value, Float_t &Bz_value, Float_t &Bphi_value ) ;
  virtual void    SetFactor (Float_t factor = 1);
  virtual void    SetRescale(Float_t factor = 1);
  virtual void    SetBDipole(Float_t m = -42.67);
  virtual void    SetRmaxDip(Float_t m =   15.3);
  virtual void    SetZminDip(Float_t m =  980.0);
  virtual void    SetZmaxDip(Float_t m = 1350.0);
  virtual void    SetLock();
  virtual EBField GetMap()     {return fMap;}
  virtual Float_t GetFactor()  {return fFactor;}
  virtual Float_t GetRescale() {return fRescale;}
  virtual Bool_t  IsLocked()   {return fLock;}
  virtual void    Print(Option_t* opt = "") const;
  void  SetStarMagFieldRotation(TGeoRotation &rot);
  void  SetStarMagFieldRotation(Double_t* rot);
  const TGeoRotation &StarMagFieldRotation() {return * &fStarMagFieldRotation;}
};

#endif
