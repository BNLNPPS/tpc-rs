#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "tpcrs/configurator.h"
#include "tpcrs/detail/mag_field.h"
#include "logger.h"


namespace tpcrs { namespace detail {


template<>
float Interpolator<float>::Linear(const float px, const float vx[2])
{
  return vx[0] + (vx[1] - vx[0]) * px;
}

template<>
float Interpolator<float>::Bilinear(const float px, const float py, const float vx[2], const float vy[2])
{
  float v_tmp[2] = { Linear(px, vx), Linear(px, vy), };
  return Linear(py, v_tmp);
}

template<>
float Interpolator<float>::Trilinear(const float px, const float py, const float pz, const float vx[2], const float vy[2], const float vz[4])
{
  float v_tmp[2] = { Bilinear(px, py, vx, vy), Bilinear(px, py, vz, vz+2) };
  return Linear(pz, v_tmp);
}


/**
 * Reads the magnetic field map from a file and scales it according to a scale
 * factor provided during instantiation. The field map accuracy = 1 G.
 *
 * This code works in kGauss, cm - but note that the field maps read from files
 * are expected to be in gauss, cm.
 */
MagField::MagField(const tpcrs::Configurator& cfg) :
  cfg_(cfg)
{
  ReadValues();
}


void MagField::ReadValues()
{
  using namespace std;

  string file_name(cfg_.S<ResponseSimulator>().mag_field_file);
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
void MagField::InterpolateField2D(const double r, const double z, double &Br, double &Bz) const
{
  using namespace std;

  int rlow = distance(begin(c1_), lower_bound(begin(c1_), end(c1_), float(r))) - 1;
  int zlow = distance(begin(c2_), lower_bound(begin(c2_), end(c2_), float(z))) - 1;

  float scale = cfg_.S<ResponseSimulator>().mag_field_scale;

  int nR = c1_.size();
  int nZ = c2_.size();

  float pr = (float(r) - c1_[rlow]) / (c1_[rlow+1] - c1_[rlow]);
  float pz = (float(z) - c2_[zlow]) / (c2_[zlow+1] - c2_[zlow]);

  Br = scale * Interpolator<float>::Bilinear(pr, pz, &v1_[zlow*nR + rlow], &v1_[(zlow + 1)*nR + rlow]);
  Bz = scale * Interpolator<float>::Bilinear(pr, pz, &v2_[zlow*nR + rlow], &v2_[(zlow + 1)*nR + rlow]);
}

} }
