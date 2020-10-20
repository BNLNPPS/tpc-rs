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
  using namespace std;

  // Scale maps to work in kGauss, cm
  float fscale = 0.001 * scale_factor_;

  float save_Br[2];
  float save_Bz[2];

  int nR = c1_.size();
  int nZ = c2_.size();

  int jlow = distance(begin(c2_), lower_bound(begin(c2_), end(c2_), float(z))) - 1;
  int klow = distance(begin(c1_), lower_bound(begin(c1_), end(c1_), float(r))) - 1;

  if (jlow < 0) jlow = 0;
  if (klow < 0) klow = 0;

  if (jlow > nZ - 2) jlow = nZ - 2;
  if (klow > nR - 2) klow = nR - 2;

  for (int j = jlow; j < jlow + 2; j++) {
    save_Br[j - jlow] = Interpolate( &c1_[klow], &v1_[j*nR + klow], r );
    save_Bz[j - jlow] = Interpolate( &c1_[klow], &v2_[j*nR + klow], r );
  }

  Br_value = fscale * Interpolate( &c2_[jlow], save_Br, z );
  Bz_value = fscale * Interpolate( &c2_[jlow], save_Bz, z );
}


/**
 * Linear interpolation
 */
float MagField::Interpolate(const float xs[2], const float ys[2], const float x)
{
  return ys[0] + ( ys[1] - ys[0] ) * ( x - xs[0] ) / ( xs[1] - xs[0] ) ;
}

} }
