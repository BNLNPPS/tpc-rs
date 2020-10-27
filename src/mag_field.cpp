#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

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
  ReadValues();
}


void MagField::ReadValues()
{
  using namespace std;

  string file_name;

  if ( field_type_ == MagFieldType::kMapped ) {                    // Mapped field values
    if ( scale_factor_ > 0 ) {
      file_name   = "bfield_full_positive_2D.dat" ;
    }
    else {
      file_name   = "bfield_full_negative_2D.dat" ;
      // The values read from file already reflect the sign
      scale_factor_ = abs(scale_factor_);
    }
  }
  else if ( field_type_ == MagFieldType::kConstant ) {           // Constant field values
    file_name = "const_full_positive_2D.dat" ;
  }

  file_name = cfg_.Locate(file_name);
  ifstream file(file_name);

  if (file)
    cout << "mag field file opened: " << file_name << '\n';
  else
    cerr << "mag field file not available: " << file_name << '\n';

  // A wrapper to read Nan and inf values from a string
  struct istringstream_ : istringstream
  {
    istringstream_(string s) : istringstream(s) {};
    istringstream_& operator>> (double& d)
    {
      string s;
      *this >> s;
      d = stod(s);
      return *this;
    }
  };

  // Temporary containers to store sorted unique values
  set<float> s1;
  set<float> s2;
  set<float> s3;

  string line; 
  while (getline(file, line))
  {
    size_t comment_indx = line.find_first_of('#');
    auto   comment_iter = (comment_indx == string::npos ? end(line) : begin(line) + comment_indx);

    line.erase(comment_iter, end(line));
    istringstream_ iss(line);

    // Skip lines with less than 4 words
    if (distance(istream_iterator<string>(iss), istream_iterator<string>()) < 4) continue;

    // First clear the stream then rewind. The order is important
    iss.clear();
    iss.seekg(0);

    double c1, c2, v1, v2;

    iss >> c1 >> c2 >> v1 >> v2;

    s1.insert(c1);
    s2.insert(c2);
    v1_.push_back(v1);
    v2_.push_back(v2);
  }

  copy(begin(s1), end(s1), back_inserter(c1_));
  copy(begin(s2), end(s2), back_inserter(c2_));
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
