#include <string>
#include <cmath>

#include "tpcrs/configurator.h"
#include "tpcrs/detail/mag_field.h"
#include "logger.h"


namespace tpcrs { namespace detail {


/**
 * Reads the magnetic field map from a file and scales it according to a scale
 * factor provided during instantiation. The field map accuracy = 1 G.
 *
 * This code works in kGauss, cm - but note that the field maps read from files
 * are expected to be in gauss, cm.
 */
MagField::MagField(const tpcrs::Configurator& cfg, MagFieldType field_type, double scale) :
  cfg_(cfg),
  field_type_(field_type),
  scale_factor_(scale)
{
  ReadField();
}


/**
 * B field in Cartesian coordinates - 3D field
 */
void MagField::GetFieldValue3D(const float x[3], float B[3]) const
{
  float Br_value = 0;
  float Bz_value = 0;
  float Bphi_value = 0;

  float r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
  float z = x[2];
  float phi;

  if ( r != 0.0 ) {
    phi = std::atan2( x[1], x[0] );

    if (phi < 0) phi += 2 * M_PI;           // Table uses phi from 0 to 2*Pi

    InterpolateField3D(r, z, phi, Br_value, Bz_value, Bphi_value);
    B[0] = Br_value * (x[0] / r) - Bphi_value * (x[1] / r);
    B[1] = Br_value * (x[1] / r) + Bphi_value * (x[0] / r);
    B[2] = Bz_value;
  }
  else {
    phi = 0 ;
    InterpolateField3D(r, z, phi, Br_value, Bz_value, Bphi_value);
    B[0] = Br_value;
    B[1] = Bphi_value;
    B[2] = Bz_value;
  }
}


/**
 * Read the electric and magnetic field maps stored on disk
 */
void MagField::ReadField()
{
  FILE*    magfile, *b3Dfile ;
  std::string filename, filename3D ;

  if ( field_type_ == MagFieldType::kMapped ) {                  	// Mapped field values
    if ( scale_factor_ > 0 ) {
      filename   = "bfield_full_positive_2D.dat" ;
      filename3D = "bfield_full_positive_3D.dat" ;
    }
    else {
      filename   = "bfield_full_negative_2D.dat" ;
      filename3D = "bfield_full_negative_3D.dat" ;
      // The values read from file already reflect the sign
      scale_factor_ = std::abs(scale_factor_);
    }
  }
  else if ( field_type_ == MagFieldType::kConstant ) {           // Constant field values
    filename = "const_full_positive_2D.dat" ;
  }

  std::string MapLocation = cfg_.Locate(filename);
  magfile = fopen(MapLocation.c_str(), "r") ;

  if (magfile)
  {
    char cname[128] ;
    fgets(cname, sizeof(cname), magfile);    // Read comment lines at begining of file
    fgets(cname, sizeof(cname), magfile);
    fgets(cname, sizeof(cname), magfile);
    fgets(cname, sizeof(cname), magfile);
    fgets(cname, sizeof(cname), magfile);

    for ( int j = 0 ; j < nZ ; j++ ) {
      for ( int k = 0 ; k < nR ; k++ ) {
        fgets  ( cname, sizeof(cname), magfile ) ;
        sscanf ( cname, " %f %f %f %f ", &R_[k], &Z_[j], &Br[j][k], &Bz[j][k] ) ;
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
  else if ( field_type_ == MagFieldType::kConstant )             // Constant field values
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
void MagField::InterpolateField2D(double r, double z, double &Br_value, double &Bz_value) const
{
  // Scale maps to work in kGauss, cm
  float fscale = 0.001 * scale_factor_;

  static int jlow = 0, klow = 0 ;
  float save_Br[2];
  float save_Bz[2];

  Search(nZ, Z_, float(z), jlow);
  Search(nR, R_, float(r), klow);

  if (jlow < 0) jlow = 0;
  if (klow < 0) klow = 0;

  if (jlow > nZ - 2) jlow = nZ - 2;
  if (klow > nR - 2) klow = nR - 2;

  for (int j = jlow; j < jlow + 2; j++) {
    save_Br[j - jlow] = Interpolate( &R_[klow], &Br[j][klow], r )   ;
    save_Bz[j - jlow] = Interpolate( &R_[klow], &Bz[j][klow], r )   ;
  }

  Br_value = fscale * Interpolate( &Z_[jlow], save_Br, z )   ;
  Bz_value = fscale * Interpolate( &Z_[jlow], save_Bz, z )   ;
}


/**
 * Interpolate the B field map - 3D interpolation
 */
void MagField::InterpolateField3D(float r, float z, float phi,
                                  float &Br_value, float &Bz_value, float &Bphi_value) const
{
  // Scale maps to work in kGauss, cm
  float fscale = 0.001 * scale_factor_;

  static  int ilow = 0, jlow = 0, klow = 0 ;
  float save_Br[2],   saved_Br[2] ;
  float save_Bz[2],   saved_Bz[2] ;
  float save_Bphi[2], saved_Bphi[2] ;

  if (r < 0) return;

  Search(nPhi, Phi3D, phi, ilow);
  Search(nZ,   Z3D,   z,   jlow);
  Search(nR,   R3D,   r,   klow);

  if ( ilow < 0 ) ilow = 0 ;
  if ( jlow < 0 ) jlow = 0 ;
  if ( klow < 0 ) klow = 0 ;

  if ( ilow > nPhi - 2 ) ilow = nPhi - 2 ;
  if ( jlow >   nZ - 2 ) jlow =   nZ - 2 ;
  if ( klow >   nR - 2 ) klow =   nR - 2 ;

  for ( int i = ilow ; i < ilow + 2 ; i++ ) {
    for ( int j = jlow ; j < jlow + 2 ; j++ ) {
      save_Br[j - jlow]   = Interpolate( &R3D[klow], &Br3D[i][j][klow], r )   ;
      save_Bz[j - jlow]   = Interpolate( &R3D[klow], &Bz3D[i][j][klow], r )   ;
      save_Bphi[j - jlow] = Interpolate( &R3D[klow], &Bphi3D[i][j][klow], r ) ;
    }

    saved_Br[i - ilow]   = Interpolate( &Z3D[jlow], save_Br, z )   ;
    saved_Bz[i - ilow]   = Interpolate( &Z3D[jlow], save_Bz, z )   ;
    saved_Bphi[i - ilow] = Interpolate( &Z3D[jlow], save_Bphi, z ) ;
  }

  Br_value   = fscale * Interpolate( &Phi3D[ilow], saved_Br, phi )   ;
  Bz_value   = fscale * Interpolate( &Phi3D[ilow], saved_Bz, phi )   ;
  Bphi_value = fscale * Interpolate( &Phi3D[ilow], saved_Bphi, phi ) ;
}


/**
 * Linear interpolation
 */
float MagField::Interpolate(const float xs[2], const float ys[2], const float x)
{
  return ys[0] + ( ys[1] - ys[0] ) * ( x - xs[0] ) / ( xs[1] - xs[0] ) ;
}


/**
 * Search an ordered table by starting at the most recently used point
 */
void MagField::Search(const int N, const float Xarray[], const float x, int &low)
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

} }
