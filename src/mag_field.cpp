#include <string>
#include <cmath>

#include "tpcrs/configurator.h"
#include "logger.h"
#include "mag_field.h"


/**
 * A package of Bfield. Methods included to read the correct Bfield map and
 * scale it according to a scale factor provided during instantiation. The
 * Bfield map accuracy FWHM = 1 G.
 *
 * The enumerations for the manual settings are:
 * enum   EBField  { kUndefined = 0, kConstant = 1, kMapped = 2 };
 * "kConstant = 1" means you wish to work with a constant, uniform, field
 * "kMapped = 2"   means you want to read values from the measured magnet maps
 *
 * This code works in kGauss, cm - but note that the Bfield maps on disk are in
 * gauss, cm.
 */
MagField::MagField(const tpcrs::Configurator& cfg, EBField map, float rescale) :
  cfg_(cfg),
  fMap(map),
  fFactor(cfg.S<MagFactor>().ScaleFactor),
  fRescale(rescale),
  fMagFieldRotation("MagFieldRotation")
{
  ReadField();
}


/**
 * B field in Cartesian coordinates - 2D field (ie. Phi symmetric)
 */
void MagField::BField( const double x[], float B[] )
{
  float Br_value = 0, Bz_value = 0, Bphi_value = 0;

  B[0] = B[1] = B[2] = 0;

  float r   = std::sqrt(x[0]*x[0] + x[1]*x[1]);
  float phi = std::atan2(x[1], x[0]);
  float z   = x[2];

  // Table uses phi from 0 to 2*Pi
  if (phi < 0) phi += 2*M_PI ;

  float za = std::abs(z);

  // within map
  if (z >= ZList[0] && z <= ZList[nZ - 1] && r <= Radius[nR - 1])
  {
    Interpolate2DBfield( r, z, Br_value, Bz_value ) ;
    double BL[3] = {0, 0, Bz_value};

    if ( r != 0.0 )      {
      BL[0] = Br_value * (x[0] / r) ;
      BL[1] = Br_value * (x[1] / r) ;
    }

    double BG[3];
    fMagFieldRotation.LocalToMaster(BL, BG);

    for (int i = 0; i < 3; i++) B[i] = BG[i];
  }
}


/// B field in Cartesian coordinates - 3D field
void MagField::B3DField( const float x[], float B[] )
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
  fMagFieldRotation.LocalToMaster(BL, BG);

  for (int i = 0; i < 3; i++) B[i] = BG[i];
}


/**
 * Read the electric and magnetic field maps stored on disk
 */
void MagField::ReadField()
{
  FILE*    magfile, *b3Dfile ;
  std::string filename, filename3D ;

  if ( fMap == kMapped ) {                  	// Mapped field values
    if ( std::abs(fFactor) > 0.8 ) {    		// Scale from full field data
      if ( fFactor > 0 ) {
        filename   = "bfield_full_positive_2D.dat" ;
        filename3D = "bfield_full_positive_3D.dat" ;
        fRescale   = 1 ;                // Normal field
      }
      else {
        filename   = "bfield_full_negative_2D.dat" ;
        filename3D = "bfield_full_negative_3D.dat" ;
        fRescale   = -1 ;               // Reversed field
      }
    }
    else {                                // Scale from half field data
      filename   = "bfield_half_positive_2D.dat" ;
      filename3D = "bfield_half_positive_3D.dat" ;
      fRescale   = 2 ;                    // Adjust scale factor to use half field data
    }
  }
  else if ( fMap == kConstant ) {           // Constant field values
    filename = "const_full_positive_2D.dat" ;
    fRescale = 1 ;                        // Normal field
  }
  else {
    LOG_ERROR << "MagField::ReadField  No map available - you must choose a mapped field or a constant field\n";
    exit(1) ;
  }

  std::string MapLocation = cfg_.Locate(filename);
  magfile = fopen(MapLocation.c_str(), "r") ;

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
      }
    }
  } else {
    LOG_ERROR << "MagField: File " << MapLocation << " not found\n";
    exit(1);
  }

  fclose(magfile) ;

  MapLocation = cfg_.Locate(filename3D);
  b3Dfile = fopen(MapLocation.c_str(), "r") ;

  if (b3Dfile)
  {
    char cname[128] ;
    // Read comment lines at begining of file
    fgets  ( cname, sizeof(cname), b3Dfile ) ;
    fgets  ( cname, sizeof(cname), b3Dfile ) ;
    fgets  ( cname, sizeof(cname), b3Dfile ) ;
    fgets  ( cname, sizeof(cname), b3Dfile ) ;
    fgets  ( cname, sizeof(cname), b3Dfile ) ;
    fgets  ( cname, sizeof(cname), b3Dfile ) ;

    for ( int i = 0 ; i < nPhi ; i++ ) {
      for ( int j = 0 ; j < nZ ; j++ ) {
        for ( int k = 0 ; k < nR ; k++ ) {
          fgets  ( cname, sizeof(cname), b3Dfile ) ;
          sscanf ( cname, " %f %f %f %f %f %f ",
                   &R3D[k], &Z3D[j], &Phi3D[i], &Br3D[i][j][k], &Bz3D[i][j][k], &Bphi3D[i][j][k] ) ;
          Phi3D[i] *= M_PI / 180. ;   // Convert to Radians  phi = 0 to 2*Pi
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
  } else {
    LOG_ERROR << "MagField: File " << MapLocation << " not found\n";
    exit(1);
  }

  fclose(b3Dfile) ;
}


/**
 * Interpolate the B field map - 2D interpolation
 */
void MagField::Interpolate2DBfield( const float r, const float z, float &Br_value, float &Bz_value )
{
  // Scale maps to work in kGauss, cm
  float fscale = 0.001 * fFactor * fRescale;

  const  int ORDER = 1; // Linear interpolation = 1, Quadratic = 2
  static int jlow = 0, klow = 0 ;
  float save_Br[ORDER + 1];
  float save_Bz[ORDER + 1];

  Search(nZ, ZList,  z, jlow);
  Search(nR, Radius, r, klow);

  if (jlow < 0) jlow = 0;
  if (klow < 0) klow = 0;

  if (jlow + ORDER >= nZ - 1) jlow = nZ - 1 - ORDER;
  if (klow + ORDER >= nR - 1) klow = nR - 1 - ORDER;

  for (int j = jlow; j < jlow + ORDER + 1; j++) {
    save_Br[j - jlow] = Interpolate( &Radius[klow], &Br[j][klow], ORDER, r )   ;
    save_Bz[j - jlow] = Interpolate( &Radius[klow], &Bz[j][klow], ORDER, r )   ;
  }

  Br_value = fscale * Interpolate( &ZList[jlow], save_Br, ORDER, z )   ;
  Bz_value = fscale * Interpolate( &ZList[jlow], save_Bz, ORDER, z )   ;
}


/**
 * Interpolate the B field map - 3D interpolation
 */
void MagField::Interpolate3DBfield( const float r, const float z, const float phi,
                                        float &Br_value, float &Bz_value, float &Bphi_value )
{
  // Scale maps to work in kGauss, cm
  float fscale = 0.001 * fFactor * fRescale;

  const   int ORDER = 1 ;                       // Linear interpolation = 1, Quadratic = 2
  static  int ilow = 0, jlow = 0, klow = 0 ;
  float save_Br[ORDER + 1],   saved_Br[ORDER + 1] ;
  float save_Bz[ORDER + 1],   saved_Bz[ORDER + 1] ;
  float save_Bphi[ORDER + 1], saved_Bphi[ORDER + 1] ;

  if (r < 0) return;

  Search(nPhi, Phi3D, phi, ilow);
  Search(nZ,   Z3D,   z,   jlow);
  Search(nR,   R3D,   r,   klow);

  if ( ilow < 0 ) ilow = 0 ;
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


/**
 * Interpolate a 3x2 table (quadratic) or a 2x2 table (linear)
 */
float MagField::Interpolate( const float Xarray[], const float Yarray[],
                                   const int ORDER, const float x )
{
  float y;

  if (ORDER == 2) // Quadratic Interpolation = 2
  {
    y  = (x - Xarray[1]) * (x - Xarray[2]) * Yarray[0] / ( (Xarray[0] - Xarray[1]) * (Xarray[0] - Xarray[2]) ) ;
    y += (x - Xarray[2]) * (x - Xarray[0]) * Yarray[1] / ( (Xarray[1] - Xarray[2]) * (Xarray[1] - Xarray[0]) ) ;
    y += (x - Xarray[0]) * (x - Xarray[1]) * Yarray[2] / ( (Xarray[2] - Xarray[0]) * (Xarray[2] - Xarray[1]) ) ;
  }
  else            // Linear Interpolation = 1
  {
    y  = Yarray[0] + ( Yarray[1] - Yarray[0] ) * ( x - Xarray[0] ) / ( Xarray[1] - Xarray[0] ) ;
  }

  return y;
}


/**
 * Search an ordered table by starting at the most recently used point
 */
void MagField::Search(int N, const float Xarray[], float x, int &low)
{
  long high ;
  int  ascend = 0, increment = 1 ;

  if ( Xarray[N - 1] >= Xarray[0] )
    ascend = 1 ; // Ascending ordered table if true

  if ( low < 0 || low > N - 1 ) {
    low = -1;
    high = N;
  }
  else  // Ordered search phase
  {
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

  // Binary search phase
  while ( (high - low) != 1 ) {
    long middle = ( high + low ) / 2 ;

    if ( (int)( x >= Xarray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x == Xarray[N - 1] ) low = N - 2 ;
  if ( x == Xarray[0]     ) low = 0 ;
}
