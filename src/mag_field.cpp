/***********************************************************************
 * Author: Jim Thomas   11/1/2000
 *
 * Description: the STAR Magnetic Field
 ***********************************************************************/

/*!

\class StarMagField

\author Jim Thomas 10 October 2000

A package of Bfield. Methods included to read the correct Bfield map and scale it
according to a scale factor provided during instantiation.
The statement from W.Love is that the Bfield map accuracy FWHM = 1 G.
<p>

An enumerated argument provided at the time of instantiation selects
a constant magnetic field (value=1) or the measured magnetic field (value=2)
at a field setting that you select manually.  Alternatively, you can use the
database to determine the magnetic field setting but you must then provide a
a time stamp and use a different instantiation (this is usually done in the chain).

The enumerations for the manual settings are:
enum   EBField  { kUndefined = 0, kConstant = 1, kMapped = 2, kChain = 3 } ;
"kConstant = 1" means you wish to work with a constant, uniform, field.
"kMapped = 2"   means you want to read values from the measured magnet maps.
The other enumerations are undefined and reserved for future expansion.

<p>

This code works in kGauss, cm - but note that the Bfield maps on disk
are in gauss, cm.

To do:  <br>
- Finish pulling parameters out of DB rather than from #define.
- Use Magnet current rather than MafFactor
- Add a routine to distort the track if we are given a Geant Vector full of points == a track
- Add simulated B field map in the regions where the field is not mapped.

*/
#include <string>
#include <cassert>
#include <cmath>

#include "TFile.h"
#include "TError.h"
#include "TEnv.h"

#include "tpcrs/configurator.h"
#include "mag_field.h"
#include "struct_containers.h"


struct BFLD_t {
  int version;
  const char* code;
  float date; int kz; float rmaxx, zmaxx, rrm, zz1, zz2;
  float RmaxInn, ZmaxInn;
  int   nrp, nzp;
};
static const BFLD_t BFLD = {// real field
  3,         // version
  "opt1",    // code:    fit version code
  22.10,     // date: 	 fit date
  22,        // kz:   	 number of z lines
  270,       // rmaxx:	 maximum radius of extrapolated measurements
  290,       // zmaxx:	 maximum length of extrapolated measurements
  400,       // rrm:  	 maximum radius of all fields
  270,       // zz1:  	 length of measured field interpolation
  800,       // zz2:  	 max length of all fields
  264.265,   // RmaxInn: Inner field volume radius
  312.500,   // ZmaxInn: Inner field volume length
  200,       // nrp:     number of R nodes in the map
  800        // nzp:     number of Z nodes in the map
};

struct BDAT_t {
  int N;
  float Zi, Ri[20], Bzi[20], Bri[20];
};

static const int nZext = 23;
static const BDAT_t BDAT[nZext] = { // calculated STAR field
  {
    15,    // Number of field points
    0.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0,
      200.0, 225.0, 250.0, 275.0, 300.0, 375.0, 400.0
    },  // Radius
    {
      4704.6, 4768.0, 4946.6, 5148.3, 5128.6, 4927.3, 4844.9, 4830.1,
      4823.3, 4858.0, 5110.9, 3402.4, -18.1, -13.6, -25.0
    },  // Axial  field
    {
      0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0
    }
  },// Radial field == 0 @ z = 0
  //    {  0.0, 131.3, 188.4,  74.1,-148.6,-164.8, -53.2,  23.8,
  //       97.4, 213.7, 329.3,  75.3,  18.2, -44.3, -36.5 } },// Radial field
  {
    15,    // Number of field points
    270.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0,
      200.0, 225.0, 250.0, 275.0, 300.0, 375.0, 400.0
    },  // Radius
    {
      4704.6, 4768.0, 4946.6, 5148.3, 5128.6, 4927.3, 4844.9, 4830.1,
      4823.3, 4858.0, 5110.9, 3402.4, -18.1, -13.6, -25.0
    },  // Axial  field
    {
      0.0, 131.3, 188.4,  74.1, -148.6, -164.8, -53.2,  23.8,
      97.4, 213.7, 329.3,  75.3,  18.2, -44.3, -36.5
    }
  },// Radial field
  {
    15,    // Number of field points
    280.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0,
      200.0, 225.0, 250.0, 275.0, 300.0, 375.0, 400.0
    },  // Radius
    {
      4568.8, 4649.8, 4898.6, 5241.5, 5234.5, 4883.9, 4806.5, 4802.4,
      4781.7, 4771.8, 5057.8, 3504.2, -144.8, -15.0, -28.0
    },  // Axial  field
    {
      0.0, 188.6, 297.5, 151.9, -241.2, -242.5, -60.1,  19.5,
      92.3, 244.3, 541.5, 396.8,  83.6, -49.9, -40.6
    }
  },// Radial field
  {
    15,    // Number of field points
    290.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0,
      200.0, 225.0, 250.0, 275.0, 300.0, 375.0, 400.0
    },  // Radius
    {
      4383.8, 4478.4, 4801.1, 5378.7, 5431.1, 4771.4, 4765.5, 4778.2,
      4741.4, 4651.9, 4852.9, 3684.4,   6.9, -16.9, -31.6
    },  // Axial  field
    {
      0.0, 260.2, 456.5, 312.9, -414.5, -349.8, -51.7,  14.4,
      74.7, 234.0, 858.0, 726.3, 355.0, -56.5, -45.0
    }
  },// Radial field
  {
    15,    // Number of field points
    300.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0,
      200.0, 225.0, 250.0, 275.0, 300.0, 375.0, 400.0
    },  // Radius
    {
      4142.1, 4240.8, 4614.8, 5546.4, 5829.1, 4450.0, 4737.6, 4761.4,
      4711.3, 4534.1, 4231.0, 4067.5, -880.0, -19.3, -36.2
    },  // Axial  field
    {
      0.0, 341.1, 669.5, 661.0, -766.7, -480.9, -24.5,   8.8,
      43.5, 149.9, 1333.6, 999.3,  53.6, -64.2, -49.8
    }
  },// Radial field
  {
    15,    // Number of field points
    310.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0,
      200.0, 225.0, 250.0, 275.0, 300.0, 375.0, 400.0
    },  // Radius
    {
      3842.7, 3930.2, 4292.5, 5589.1, 6643.0, 3236.8, 4733.0, 4755.1,
      4699.4, 4485.0, 1931.8, 4782.0,  50.2, -22.8, -42.0
    },  // Axial  field
    {
      0.0, 421.2, 915.6, 1382.6, -1482.8, -1019.7,  1.2,  2.0,
      1.9,  -2.3, 2069.4, 791.7, 240.6, -73.6, -54.9
    }
  },// Radial field
  {
    8,    // Number of field points
    320.0,    // distance to  Z=0
    {  0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 375.0, 400.0 },  // Radius
    {3491.2, 3552.1, 3807.3, 4923.7, 7889.6, 1983.9, -28.0, -49.4 }, // Axial  field
    {  0.0, 485.7, 1133.5, 2502.8, -38.8, -174.8, -85.1, -60.0 }
  },// Radial field
  {
    7,    // Number of field points
    330.0,    // distance to  Z=0
    {  0.0,  25.0,  50.0,  75.0, 250.0, 375.0, 400.0 },  // Radius
    {3105.3, 3127.0, 3200.4, 3268.9,  -3.5, -36.6, -59.0 }, // Axial  field
    {  0.0, 521.1, 1246.1, 3029.5, 9199.2, -99.4, -64.5 }
  },// Radial field
  {
    6,    // Number of field points
    340.0,    // distance to  Z=0
    {  0.0,  25.0,  50.0,  75.0, 375.0, 400.0 },  // Radius
    {2706.4, 2686.8, 2574.5, 1826.7, -51.8, -71.0 }, // Axial  field
    {  0.0, 520.6, 1218.1, 2485.3, -116.9, -67.3 }
  },// Radial field
  {
    6,    // Number of field points
    350.0,    // distance to  Z=0
    {  0.0,  25.0,  50.0,  75.0, 375.0, 400.0 },  // Radius
    {2317.7, 2264.6, 2026.3, 1142.6, -80.8, -85.1 }, // Axial  field
    {  0.0, 487.6, 1082.3, 1787.2, -133.8, -67.0 }
  },// Radial field
  {
    8,    // Number of field points
    360.0,    // distance to  Z=0
    {  0.0,  25.0,  50.0,  75.0, 100.0, 250.0, 375.0, 400.0 },  // Radius
    {1958.5, 1885.6, 1595.6, 829.2, -563.7, 4895.8, -127.6, -99.8 }, // Axial  field
    {  0.0, 432.4, 901.7, 1265.8, 788.0, 9507.4, -134.0, -62.2 }
  },// Radial field
  {
    17,    // Number of field points
    370.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0, 200.0,
      225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0
    },  // Radius
    {
      1637.8, 1562.2, 1276.7, 678.9,  15.7, 251.6, 384.9, 503.7, 683.3,
      1087.1, 1868.1, -1320.5, -593.9, -391.5, -345.9, -168.2, -112.9
    },  // Axial field
    {
      0.0, 367.3, 720.6, 900.1, 421.6,  60.4,  37.1,  44.5,  79.7,
      229.6, 2339.4, 654.6, 114.6,  35.9, -30.0, -101.8, -52.4
    }
  },// Radial field
  {
    17,    // Number of field points
    380.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0, 200.0,
      225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0
    },  // Radius
    {
      1373.5, 1296.6, 1045.5, 603.2, 221.4, 278.6, 382.7, 488.2, 638.7,
      892.4, 708.6, -709.9, -515.0, -364.7, -293.1, -181.5, -122.1
    },  // Axial  field
    {
      0.0, 302.3, 563.3, 650.3, 369.7, 120.0,  79.6,  96.2, 169.1,
      430.1, 1454.7, 860.7, 228.6,  77.5, -10.8, -60.2, -39.4
    }
  },// Radial field
  {
    17,    // Number of field points
    390.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0, 200.0,
      225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0
    },  // Radius
    {
      1151.2, 1083.6, 877.2, 557.6, 308.9, 305.2, 377.6, 463.3, 573.3,
      684.5, 377.5, -376.2, -415.2, -326.0, -258.2, -179.8, -126.9
    },  // Axial  field
    {
      0.0, 243.7, 437.7, 486.1, 319.9, 155.1, 115.2, 139.4, 232.4,
      494.6, 1019.1, 751.4, 289.6, 112.2,  19.4, -26.7, -25.0
    }
  },// Radial field
  {
    17,    // Number of field points
    400.0,    // distance to  Z=0
    {
      0.0,  25.0,  50.0,  75.0, 100.0, 125.0, 150.0, 175.0, 200.0,
      225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 375.0, 400.0
    },  // Radius
    {
      971.6, 914.8, 751.6, 520.8, 348.0, 323.7, 369.1, 432.0, 500.9,
      520.4, 251.2, -214.3, -320.8, -282.3, -230.2, -171.7, -127.7
    },  // Axial  field
    {
      0.0, 194.5, 341.1, 375.8, 277.5, 171.7, 142.1, 172.1, 269.4,
      486.6, 769.0, 624.1, 308.0, 137.2,  44.9,  -1.4, -11.2
    }
  },// Radial field
  {
    9,    // Number of field points
    450.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    {481.5, 423.2, 325.1, 283.8, 242.6,  88.2, -85.2, -123.1, -103.9 }, // Axial
    {  0.0, 119.3, 157.5, 174.6, 248.4, 314.8, 220.8,  95.4,  32.3 }
  },// Radial
  {
    9,    // Number of field points
    500.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    {291.3, 273.2, 234.4, 192.6, 136.7,  53.6, -25.6, -64.6, -72.5 },  // Axial
    {  0.0,  60.7, 103.2, 135.9, 168.1, 177.4, 140.6,  84.7,  41.9 }
  },// Radial
  {
    9,    // Number of field points
    550.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    {190.4, 181.9, 159.6, 127.7,  86.0,  37.1,  -7.3, -35.9, -48.9 },  // Axial
    {  0.0,  37.8,  69.7,  94.3, 110.3, 110.8,  92.6,  64.2,  37.3 }
  },// Radial
  {
    9,    // Number of field points
    600.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    {126.9, 121.7, 107.1,  84.9,  56.8,  26.4,  -1.2, -21.2, -32.9 },  // Axial
    {  0.0,  25.0,  46.9,  63.4,  72.6,  72.2,  62.3,  46.3,  29.2 }
  },// Radial
  {
    9,    // Number of field points
    650.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    { 84.8,  81.4,  71.7,  56.8,  38.3,  18.7,   0.8, -13.1, -22.2 },  // Axial
    {  0.0,  16.7,  31.4,  42.4,  48.2,  48.1,  42.4,  32.8,  21.6 }
  },// Radial
  {
    9,    // Number of field points
    700.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    { 56.7,  54.4,  47.9,  38.0,  25.9,  13.1,   1.3,  -8.3, -15.0 },  // Axial
    {  0.0,  11.2,  21.1,  28.4,  32.4,  32.5,  29.1,  23.0,  15.5 }
  },// Radial
  {
    9,    // Number of field points
    750.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    { 37.7,  36.2,  31.9,  25.4,  17.5,   9.1,   1.2,  -5.4, -10.1 },  // Axial
    {  0.0,   7.6,  14.3,  19.3,  22.0,  22.2,  20.1,  16.1,  11.0 }
  },// Radial
  {
    9,    // Number of field points
    800.0,    // distance to  Z=0
    {  0.0,  50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0 },  // Radius
    { 24.8,  23.8,  21.0,  16.8,  11.6,   6.1,   0.9,  -3.5,  -6.7 },  // Axial
    {  0.0,   5.2,   9.8,  13.3,  15.2,  15.4,  14.1,  11.4,   7.9 }
  },// Radial
};


bool StarMagField::mConstBz = false;



StarMagField::StarMagField ( EBField map, float factor,
                             bool lock, float rescale,
                             float BDipole, float RmaxDip,
                             float ZminDip, float ZmaxDip) :
  fBzdZCorrection(0),
  fBrdZCorrection(0),
  fMap(map),
  fFactor(factor),   fRescale(rescale),
  fBDipole(BDipole), fRmaxDip(RmaxDip),
  fZminDip(ZminDip), fZmaxDip(ZmaxDip),
  fLock(lock)
{
  if (fMap == kUndefined) {
    printf("StarMagField is instantiated with predefined factor %f and map %i\n", fFactor, fMap);
  }
  else {
    if (fLock) printf("StarMagField is locked, no modification from DB will be accepted\n");
  }

  ReadField() ;                       // Read the Magnetic
  fStarMagFieldRotation = TGeoRotation("StarMagFieldRotation");
}


/// B field in Cartesian coordinates - 2D field (ie. Phi symmetric)
void StarMagField::BField( const double x[], double B[] )
{
  float xx[3] = {(float) x[0], (float) x[1], (float) x[2]};
  float bb[3];
  BField(xx, bb);
  B[0] = bb[0]; B[1] = bb[1]; B[2] = bb[2];
}


void StarMagField::BField( const float x[], float B[] )

{


  float r, z, Br_value, Bz_value ;
  float phi, Bphi_value, phi1;
  Bphi_value = 0;
  Br_value =  Bz_value = 0;
  B[0] = B[1] = B[2] = 0;
  z  = x[2] ;
  r  = std::sqrt( x[0] * x[0] + x[1] * x[1] ) ;
  phi = std::atan2( x[1], x[0] ) ;

  if ( phi < 0 ) phi += 2 * M_PI ;           // Table uses phi from 0 to 2*Pi



  if ( mConstBz ) {
    B[0] = B[1] = B[2] = 0.;

    if ( abs(z) < 380.0 && r < 300.0 ) B[2] = +5.0;

    return;
  }



  float za = std::abs(z);

  if (za > fZminDip && za < fZmaxDip && r < fRmaxDip) {//     beam Dipole
    B[1] = std::copysign(fBDipole, z);
    B[2] = std::abs(B[1] / 1000.);
    return;
  }

  if (z >= ZList[0] && z <= ZList[nZ - 1] && r <= Radius[nR - 1]) { // within Map
    Interpolate2DBfield( r, z, Br_value, Bz_value ) ;
    double BL[3] = {0, 0, Bz_value};

    if ( r != 0.0 )      {
      BL[0] = Br_value * (x[0] / r) ;
      BL[1] = Br_value * (x[1] / r) ;
    }

    double BG[3];
    fStarMagFieldRotation.LocalToMaster(BL, BG);

    for (int i = 0; i < 3; i++) B[i] = BG[i];

    return;
  }

  //   //added by Lijuan within the steel


  if (za <= 342.20  && r >= 303.29 && r <= 364.25) { // within Map

    phi1 = phi * 180.0 / M_PI;

    if (phi1 > 12) phi1 = phi1 - int(phi1 / 12) * 12;

    Interpolate3DBSteelfield( r, za, phi1, Br_value, Bz_value, Bphi_value ) ;
    B[0] = Br_value * (x[0] / r) - Bphi_value * (x[1] / r) ;
    B[1] = Br_value * (x[1] / r) + Bphi_value * (x[0] / r) ;
    B[2] = Bz_value ;

    if (z < 0) {
      B[0] = -B[0];
      B[1] = -B[1];
    }

    return;
  }

  //end added by Lijuan within the steel


  Interpolate2ExtDBfield( r, z, Br_value, Bz_value ) ;

  if (za <= BFLD.zmaxx && r <= BFLD.rmaxx) {
    static const float zero = 0;
    static const float one = 1;
    float wz = (za - ZList[nZ - 1] ) / (BFLD.zmaxx - ZList[nZ - 1]);
    float wr = (r  - Radius[nR - 1]) / (BFLD.rmaxx - Radius[nR - 1]);
    float w  = std::min(std::max(zero, std::max(wz, wr)), one);
    float rm = std::min(r, Radius[nR - 1]);
    float zm = std::copysign(std::min(za, ZList[nZ - 1]), z);
    float BrI, BzI;
    Interpolate2DBfield( rm, zm, BrI, BzI ) ;
    Br_value = (1 - w) * BrI + w * Br_value;
    Bz_value = (1 - w) * BzI + w * Bz_value;
  }

  B[2] = Bz_value ;

  if ( r != 0.0 )      {
    B[0] = Br_value * (x[0] / r) ;
    B[1] = Br_value * (x[1] / r) ;
  }

  return;
}


/// Bfield in Cartesian coordinates - 3D field
void StarMagField::B3DField( const float x[], float B[] )
{
  float r, z, phi, Br_value, Bz_value, Bphi_value ;
  Bphi_value = 0;
  Br_value =  Bz_value = 0;
  B[0] = B[1] = B[2] = 0;
  z = x[2] ;
  r  = sqrt( x[0] * x[0] + x[1] * x[1] ) ;

  if ( r != 0.0 ) {
    phi = std::atan2( x[1], x[0] ) ;

    if ( phi < 0 ) phi += 2 * M_PI ;           // Table uses phi from 0 to 2*Pi

    Interpolate3DBfield( r, z, phi, Br_value, Bz_value, Bphi_value ) ;
    B[0] = Br_value * (x[0] / r) - Bphi_value * (x[1] / r) ;
    B[1] = Br_value * (x[1] / r) + Bphi_value * (x[0] / r) ;
    B[2] = Bz_value ;
  }
  else {
    phi = 0 ;
    Interpolate3DBfield( r, z, phi, Br_value, Bz_value, Bphi_value ) ;
    B[0] = Br_value ;
    B[1] = Bphi_value ;
    B[2] = Bz_value ;
  }

  double BL[3] = {B[0], B[1], B[2]};
  double BG[3];
  fStarMagFieldRotation.LocalToMaster(BL, BG);

  for (int i = 0; i < 3; i++) B[i] = BG[i];

  return ;

}


/// Read the electric and magnetic field maps stored on disk

void StarMagField::ReadField( )

{
  FILE*    magfile, *b3Dfile ;
  std::string comment, filename, filename3D ;
  std::string MapLocation ;

  if (gEnv->GetValue("NewTpcAlignment", 0) != 0) {
    TFile pFile(tpcrs::Configurator::Locate("StarFieldZ.root").c_str());
    TH2F* Br0 = (TH2F*) pFile.Get("Br0");
    TH2F* Bz0 = (TH2F*) pFile.Get("Bz0");

    if (Br0 && Bz0) {
      TH2F* Br5cm = (TH2F*) pFile.Get("Br5cm");
      TH2F* Bz5cm = (TH2F*) pFile.Get("Bz5cm");
      assert(Br5cm && Bz5cm);
      TH2F* Br10cm = (TH2F*) pFile.Get("Br10cm");
      TH2F* Bz10cm = (TH2F*) pFile.Get("Bz10cm");
      assert(Br10cm && Bz10cm);
      fBzdZCorrection = new TH2F(*Bz5cm); fBzdZCorrection->SetDirectory(0);
      fBzdZCorrection->Scale(0.5);
      fBzdZCorrection->Add(Bz10cm, 0.5);
      fBzdZCorrection->Add(Bz0, -1.0);
      fBrdZCorrection = new TH2F(*Br5cm); fBrdZCorrection->SetDirectory(0);
      fBrdZCorrection->Scale(0.5);
      fBrdZCorrection->Add(Br10cm, 0.5);
      fBrdZCorrection->Add(Br0, -1.0);
      Warning("StarMagField::ReadField", "Use effective PMT box dZ = 7.5 cm");
    }
  }

  if ( fMap == kMapped ) {                  	// Mapped field values
    if ( std::abs(fFactor) > 0.8 ) {    		// Scale from full field data
      if ( fFactor > 0 ) {
        filename   = "bfield_full_positive_2D.dat" ;
        filename3D = "bfield_full_positive_3D.dat" ;
        comment    = "Measured Full Field" ;
        fRescale   = 1 ;                // Normal field
      }
      else {
        filename   = "bfield_full_negative_2D.dat" ;
        filename3D = "bfield_full_negative_3D.dat" ;
        comment    = "Measured Full Field Reversed" ;
        fRescale   = -1 ;               // Reversed field
      }
    }
    else {                                // Scale from half field data
      filename   = "bfield_half_positive_2D.dat" ;
      filename3D = "bfield_half_positive_3D.dat" ;
      comment    = "Measured Half Field" ;
      fRescale   = 2 ;                    // Adjust scale factor to use half field data
    }
  }
  else if ( fMap == kConstant ) {           // Constant field values
    filename = "const_full_positive_2D.dat" ;
    comment  = "Constant Full Field" ;
    fRescale = 1 ;                        // Normal field
  }
  else {
    fprintf(stderr, "StarMagField::ReadField  No map available - you must choose a mapped field or a constant field\n");
    exit(1) ;
  }

  printf("StarMagField::ReadField  Reading  Magnetic Field  %s,  Scale factor = %f \n", comment.c_str(), fFactor);
  printf("StarMagField::ReadField  Filename is %s, Adjusted Scale factor = %f \n", filename.c_str(), fFactor * fRescale);

  MapLocation = tpcrs::Configurator::Locate(filename);
  magfile = fopen(MapLocation.c_str(), "r") ;
  printf("StarMagField::ReadField  Reading  2D Magnetic Field file: %s \n", filename.c_str());

  if (magfile)

  {
    char cname[128] ;
    fgets  ( cname, sizeof(cname), magfile ) ;     // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), magfile ) ;
    fgets  ( cname, sizeof(cname), magfile ) ;
    fgets  ( cname, sizeof(cname), magfile ) ;
    fgets  ( cname, sizeof(cname), magfile ) ;

    for ( int j = 0 ; j < nZ ; j++ ) {
      for ( int k = 0 ; k < nR ; k++ ) {
        fgets  ( cname, sizeof(cname), magfile ) ;
        sscanf ( cname, " %f %f %f %f ", &Radius[k], &ZList[j], &Br[j][k], &Bz[j][k] ) ;

        if (fBzdZCorrection && fBrdZCorrection) {
          Br[j][k] += fFactor * fBrdZCorrection->Interpolate(ZList[j], Radius[k]);
          Bz[j][k] += fFactor * fBzdZCorrection->Interpolate(ZList[j], Radius[k]);
        }
      }
    }
  }

  else

  {
    fprintf(stderr, "StarMagField::ReadField  File %s not found !\n", MapLocation.c_str());
    exit(1);
  }

  fclose(magfile) ;

  MapLocation = tpcrs::Configurator::Locate(filename3D);
  b3Dfile = fopen(MapLocation.c_str(), "r") ;
  printf("StarMagField::ReadField  Reading 3D Magnetic Field file: %s \n", filename3D.c_str());

  if (b3Dfile)

  {
    char cname[128] ;
    fgets  ( cname, sizeof(cname), b3Dfile ) ;     // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), b3Dfile ) ;     // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), b3Dfile ) ;     // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), b3Dfile ) ;     // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), b3Dfile ) ;     // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), b3Dfile ) ;     // Read comment lines at begining of file

    for ( int i = 0 ; i < nPhi ; i++ ) {
      for ( int j = 0 ; j < nZ ; j++ ) {
        for ( int k = 0 ; k < nR ; k++ ) {
          fgets  ( cname, sizeof(cname), b3Dfile ) ;
          sscanf ( cname, " %f %f %f %f %f %f ",
                   &R3D[k], &Z3D[j], &Phi3D[i], &Br3D[i][j][k], &Bz3D[i][j][k], &Bphi3D[i][j][k] ) ;
          Phi3D[i] *= M_PI / 180. ;   // Convert to Radians  phi = 0 to 2*Pi

          if (fBzdZCorrection && fBrdZCorrection) {
            Br3D[i][j][k] += fFactor * fBrdZCorrection->Interpolate(Z3D[j], R3D[k]);
            Bz3D[i][j][k] += fFactor * fBzdZCorrection->Interpolate(Z3D[j], R3D[k]);
          }
        }
      }
    }
  }

  else if ( fMap == kConstant )             // Constant field values

  {
    for ( int i = 0 ; i < nPhi ; i++ ) {
      for ( int j = 0 ; j < nZ ; j++ ) {
        for ( int k = 0 ; k < nR ; k++ ) {
          Br3D[i][j][k] = Br[j][k] ;
          Bz3D[i][j][k] = Bz[j][k] ;
          Bphi3D[i][j][k] = 0 ;
        }
      }
    }
  }

  else

  {
    fprintf(stderr, "StarMagField::ReadField  File %s not found !\n", MapLocation.c_str());
    exit(1);
  }

  fclose(b3Dfile) ;
  MapLocation = tpcrs::Configurator::Locate("steel_magfieldmap.dat");
  magfile = fopen(MapLocation.c_str(), "r") ;

  if (magfile) {
    printf("StarMagField::ReadField  Reading  3D Magnetic Field file: %s \n", filename.c_str());
    char cname[128] ;

    for (;;) {
      fgets  ( cname, sizeof(cname), magfile ) ;     // Read comment lines at begining of file

      if (cname[0] == '#') continue;

      break;
    }

    for ( int i = 0 ; i < nPhiSteel ; i++ ) {

      for ( int k = 0 ; k < nRSteel ; k++ ) {
        for ( int j = 0 ; j < nZSteel ; j++ ) {


          fgets  ( cname, sizeof(cname), magfile ) ;
          sscanf ( cname, " %f %f %f %f %f %f ",
                   &R3DSteel[k], &Z3DSteel[j], &Phi3DSteel[i], &Bx3DSteel[i][j][k], &Bz3DSteel[i][j][k], &By3DSteel[i][j][k] ) ;

          //added by Lijuan
          Br3DSteel[i][j][k] = std::cos(Phi3DSteel[i] * M_PI / 180.) * Bx3DSteel[i][j][k] + std::sin(Phi3DSteel[i] * M_PI / 180.) * By3DSteel[i][j][k];

          Bphi3DSteel[i][j][k] = 0 - std::sin(Phi3DSteel[i] * M_PI / 180.) * Bx3DSteel[i][j][k] + std::cos(Phi3DSteel[i] * M_PI / 180.) * By3DSteel[i][j][k];
        }
      }
    }

    fclose(magfile);
  }

  return ;

}





/// Interpolate the B field map - 2D interpolation

void StarMagField::Interpolate2DBfield( const float r, const float z, float &Br_value, float &Bz_value )

{

  float fscale ;

  fscale = 0.001 * fFactor * fRescale ;           // Scale STAR maps to work in kGauss, cm


  const   int ORDER = 1  ;                      // Linear interpolation = 1, Quadratic = 2
  static  int jlow = 0, klow = 0 ;
  float save_Br[ORDER + 1] ;
  float save_Bz[ORDER + 1] ;

  Search ( nZ, ZList,  z, jlow ) ;
  Search ( nR, Radius, r, klow ) ;

  if ( jlow < 0 ) jlow = 0 ;   // artifact of Root's binsearch, returns -1 if out of range

  if ( klow < 0 ) klow = 0 ;

  if ( jlow + ORDER  >=    nZ - 1 ) jlow =   nZ - 1 - ORDER ;

  if ( klow + ORDER  >=    nR - 1 ) klow =   nR - 1 - ORDER ;

  for ( int j = jlow ; j < jlow + ORDER + 1 ; j++ ) {
    save_Br[j - jlow]   = Interpolate( &Radius[klow], &Br[j][klow], ORDER, r )   ;
    save_Bz[j - jlow]   = Interpolate( &Radius[klow], &Bz[j][klow], ORDER, r )   ;
  }

  Br_value  = fscale * Interpolate( &ZList[jlow], save_Br, ORDER, z )   ;
  Bz_value  = fscale * Interpolate( &ZList[jlow], save_Bz, ORDER, z )   ;

}


void StarMagField::Interpolate2ExtDBfield( const float r, const float z, float &Br_value, float &Bz_value )
{
  static float ZExtList[nZext];
  static bool  first = true;

  if (first) {
    for (int j = 0; j < nZext; j++) ZExtList[j] = BDAT[j].Zi;

    first = false;
  }

  float za = std::abs(z);

  if (za > BFLD.zz2 || r > BFLD.rrm) return;

  if (za < ZList[nZ - 1] && r < Radius[nR - 1]) return;

  //added by Lijuan
  if (za <= 342.20  && r >= 303.29 && r <= 363.29) return;

  //end added by Lijuan


  float fscale  = 0.001 * fFactor; // Scale STAR maps to work in kGauss, cm. Table only for Full Field, no Rescale !

  const   int ORDER = 1  ;                      // Linear interpolation = 1, Quadratic = 2
  static  int jlow = 0, klow = 0 ;
  float save_Br[ORDER + 1] ;
  float save_Bz[ORDER + 1] ;
  Search ( nZext, ZExtList,  za, jlow ) ;

  if ( jlow < 0 ) jlow = 0 ;   // artifact of Root's binsearch, returns -1 if out of range

  if ( jlow + ORDER  >=    nZext - 1 ) jlow =   nZext - 1 - ORDER ;

  for ( int j = jlow ; j < jlow + ORDER + 1 ; j++ ) {
    int N = BDAT[j].N;
    Search ( N, (float*) (&BDAT[j].Ri[0]), r, klow ) ;

    if ( klow < 0 ) klow = 0 ;

    if ( klow + ORDER  >=    BDAT[j].N - 1 ) klow =   BDAT[j].N - 1 - ORDER ;

    save_Br[j - jlow]   = Interpolate( &BDAT[j].Ri[klow], &BDAT[j].Bri[klow], ORDER, r )   ;
    save_Bz[j - jlow]   = Interpolate( &BDAT[j].Ri[klow], &BDAT[j].Bzi[klow], ORDER, r )   ;
  }

  Br_value  = fscale * Interpolate( &ZExtList[jlow], save_Br, ORDER, za )   ;
  Bz_value  = fscale * Interpolate( &ZExtList[jlow], save_Bz, ORDER, za )   ;

  if (z < 0) Br_value  = - Br_value;
}

/// Interpolate the B field map - 3D interpolation

void StarMagField::Interpolate3DBfield( const float r, const float z, const float phi,
                                        float &Br_value, float &Bz_value, float &Bphi_value )
{

  float fscale ;

  fscale = 0.001 * fFactor * fRescale ;           // Scale STAR maps to work in kGauss, cm

  const   int ORDER = 1 ;                       // Linear interpolation = 1, Quadratic = 2
  static  int ilow = 0, jlow = 0, klow = 0 ;
  float save_Br[ORDER + 1],   saved_Br[ORDER + 1] ;
  float save_Bz[ORDER + 1],   saved_Bz[ORDER + 1] ;
  float save_Bphi[ORDER + 1], saved_Bphi[ORDER + 1] ;



  //LOG_INFO<<"r===  "<<r<<"  z===  "<<z<<"  phi===  "<<phi<<'\n';
  if (r < 0) return;

  Search( nPhi, Phi3D, phi, ilow ) ;
  Search( nZ,   Z3D,   z,   jlow ) ;
  Search( nR,   R3D,   r,   klow ) ;

  if ( ilow < 0 ) ilow = 0 ;   // artifact of Root's binsearch, returns -1 if out of range

  if ( jlow < 0 ) jlow = 0 ;

  if ( klow < 0 ) klow = 0 ;

  if ( ilow + ORDER  >=  nPhi - 1 ) ilow = nPhi - 1 - ORDER ;

  if ( jlow + ORDER  >=    nZ - 1 ) jlow =   nZ - 1 - ORDER ;

  if ( klow + ORDER  >=    nR - 1 ) klow =   nR - 1 - ORDER ;

  for ( int i = ilow ; i < ilow + ORDER + 1 ; i++ ) {
    for ( int j = jlow ; j < jlow + ORDER + 1 ; j++ ) {
      save_Br[j - jlow]   = Interpolate( &R3D[klow], &Br3D[i][j][klow], ORDER, r )   ;
      save_Bz[j - jlow]   = Interpolate( &R3D[klow], &Bz3D[i][j][klow], ORDER, r )   ;
      save_Bphi[j - jlow] = Interpolate( &R3D[klow], &Bphi3D[i][j][klow], ORDER, r ) ;
    }

    saved_Br[i - ilow]   = Interpolate( &Z3D[jlow], save_Br, ORDER, z )   ;
    saved_Bz[i - ilow]   = Interpolate( &Z3D[jlow], save_Bz, ORDER, z )   ;
    saved_Bphi[i - ilow] = Interpolate( &Z3D[jlow], save_Bphi, ORDER, z ) ;
  }

  Br_value   = fscale * Interpolate( &Phi3D[ilow], saved_Br, ORDER, phi )   ;
  Bz_value   = fscale * Interpolate( &Phi3D[ilow], saved_Bz, ORDER, phi )   ;
  Bphi_value = fscale * Interpolate( &Phi3D[ilow], saved_Bphi, ORDER, phi ) ;

}








//added by Lijuan for the magnetic field in steel.


/// Interpolate the B field map - 3D interpolation

void StarMagField::Interpolate3DBSteelfield( const float r, const float z, const float phi,
    float &Br_value, float &Bz_value, float &Bphi_value )
{

  float fscale ;

  //This is different from the usual bfield map, changed by Lijuan

  //   fscale = 0.001*fFactor*fRescale ;               // Scale STAR maps to work in kGauss, cm
  fscale = 0.001 * fFactor;             // Scale STAR maps to work in kGauss, cm

  const   int ORDER = 1 ;                       // Linear interpolation = 1, Quadratic = 2
  static  int ilow = 0, jlow = 0, klow = 0 ;
  float save_Br[ORDER + 1],   saved_Br[ORDER + 1] ;
  float save_Bz[ORDER + 1],   saved_Bz[ORDER + 1] ;
  float save_Bphi[ORDER + 1], saved_Bphi[ORDER + 1] ;
  //  phi=phi+1;

  Search( nPhiSteel, Phi3DSteel, phi, ilow ) ;
  Search( nZSteel,   Z3DSteel,   z,   jlow ) ;
  Search( nRSteel,   R3DSteel,   r,   klow ) ;

  if ( ilow < 0 ) ilow = 0 ;   // artifact of Root's binsearch, returns -1 if out of range

  if ( jlow < 0 ) jlow = 0 ;

  if ( klow < 0 ) klow = 0 ;

  if ( ilow + ORDER  >=  nPhiSteel - 1 ) ilow = nPhiSteel - 1 - ORDER ;

  if ( jlow + ORDER  >=    nZSteel - 1 ) jlow =   nZSteel - 1 - ORDER ;

  if ( klow + ORDER  >=    nRSteel - 1 ) klow =   nRSteel - 1 - ORDER ;

  for ( int i = ilow ; i < ilow + ORDER + 1 ; i++ ) {
    for ( int j = jlow ; j < jlow + ORDER + 1 ; j++ ) {
      save_Br[j - jlow]   = Interpolate( &R3DSteel[klow], &Br3DSteel[i][j][klow], ORDER, r )   ;
      save_Bz[j - jlow]   = Interpolate( &R3DSteel[klow], &Bz3DSteel[i][j][klow], ORDER, r )   ;
      save_Bphi[j - jlow] = Interpolate( &R3DSteel[klow], &Bphi3DSteel[i][j][klow], ORDER, r ) ;
    }

    saved_Br[i - ilow]   = Interpolate( &Z3DSteel[jlow], save_Br, ORDER, z )   ;
    saved_Bz[i - ilow]   = Interpolate( &Z3DSteel[jlow], save_Bz, ORDER, z )   ;
    saved_Bphi[i - ilow] = Interpolate( &Z3DSteel[jlow], save_Bphi, ORDER, z ) ;
  }

  Br_value   = fscale * Interpolate( &Phi3DSteel[ilow], saved_Br, ORDER, phi )   ;
  Bz_value   = fscale * Interpolate( &Phi3DSteel[ilow], saved_Bz, ORDER, phi )   ;
  Bphi_value = fscale * Interpolate( &Phi3DSteel[ilow], saved_Bphi, ORDER, phi ) ;

}


/// Interpolate a 3x2 table (quadratic) or a 2x2 table (linear)
float StarMagField::Interpolate( const float Xarray[], const float Yarray[],
                                   const int ORDER, const float x )

{

  float y ;


  if ( ORDER == 2 )                // Quadratic Interpolation = 2

  {
    y  = (x - Xarray[1]) * (x - Xarray[2]) * Yarray[0] / ( (Xarray[0] - Xarray[1]) * (Xarray[0] - Xarray[2]) ) ;
    y += (x - Xarray[2]) * (x - Xarray[0]) * Yarray[1] / ( (Xarray[1] - Xarray[2]) * (Xarray[1] - Xarray[0]) ) ;
    y += (x - Xarray[0]) * (x - Xarray[1]) * Yarray[2] / ( (Xarray[2] - Xarray[0]) * (Xarray[2] - Xarray[1]) ) ;

  }

  else                             // Linear Interpolation = 1

  {
    y  = Yarray[0] + ( Yarray[1] - Yarray[0] ) * ( x - Xarray[0] ) / ( Xarray[1] - Xarray[0] ) ;
  }

  return (y) ;

}





/// Search an ordered table by starting at the most recently used point

void StarMagField::Search( int N, const float Xarray[], float x, int &low )

{
  assert(! std::isnan(x));
  long middle, high ;
  int  ascend = 0, increment = 1 ;

  if ( Xarray[N - 1] >= Xarray[0] ) ascend = 1 ; // Ascending ordered table if true

  if ( low < 0 || low > N - 1 ) { low = -1 ; high = N ; }

  else {                                          // Ordered Search phase
    if ( (int)( x >= Xarray[low] ) == ascend ) {
      if ( low == N - 1 ) return ;

      high = low + 1 ;

      while ( (int)( x >= Xarray[high] ) == ascend ) {
        low = high ;
        increment *= 2 ;
        high = low + increment ;

        if ( high > N - 1 )  {  high = N ; break ;  }
      }
    }
    else {
      if ( low == 0 )  {  low = -1 ;  return ;  }

      high = low - 1 ;

      while ( (int)( x < Xarray[low] ) == ascend ) {
        high = low ;
        increment *= 2 ;

        if ( increment >= high )  {  low = -1 ;  break ;  }
        else  low = high - increment ;
      }
    }
  }

  while ( (high - low) != 1 ) {                  // Binary Search Phase
    middle = ( high + low ) / 2 ;

    if ( (int)( x >= Xarray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x == Xarray[N - 1] ) low = N - 2 ;

  if ( x == Xarray[0]   ) low = 0 ;

}
