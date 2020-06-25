#ifndef TPCRS_TPCRS_H_
#define TPCRS_TPCRS_H_

#include <vector>

#include "tpcrs/configurator.h"

namespace tpcrs {

/**
 * Associates a packed hardware index with the data.
 *
 *                sector   row    pad   timebin
 * projected max      24    72    182       512
 * bits                5 +   7 +   10 +      10 = 32 = 4 bytes = unsigned int
 * max value          32   128   1024      1024
 */
struct DigiChannel
{
  unsigned int sector : 5, row : 7, pad : 10, timebin : 10;
  short adc;
  short idt;
};


struct SimuHit
{
  int np;         /// Number of primary electrons
  double de;      /// Energy deposited at hit
  double ds;      /// path length within pad row
  float adc;
  float pad;
  float timebin;
};


struct SimulatedCharge
{
  unsigned int sector : 5, row : 7, pad : 10, timebin : 10;
  double q;
};


struct GeantHit
{
  int track_id;    /// Parent track
  int particle_id; /// GEANT particle id
  int volume_id;   /// Volume id packed as SSRR, SS = sector, RR = pad row
  double x[3];     /// Hit center
  double p[3];     /// Local momentum
  double de;       /// energy deposited at hit
  double ds;       /// path length within pad row
  double s;        /// track length up to this hit. Used in hit ordering
  double tof;      /// Time of flight including the GEANT vertex production time
  float lgam;      /// Deprecated. ALOG10(GEKin/AMass)
  SimuHit digi;    /// Deprecated.
};


struct TPC
{
  enum Half { first, second };
};


/**
 * Expects a type similar to GeneratedHit
 */
template<typename T>
inline bool operator< (const T& lhs, const T& rhs)
{
  // sectors
  if ((lhs.volume_id % 100000) / 100 != (rhs.volume_id % 100000) / 100)
    return (lhs.volume_id % 100000) / 100 < (rhs.volume_id % 100000) / 100;

  // track id
  if (lhs.track_id != rhs.track_id)
    return lhs.track_id < rhs.track_id;

  // track length
  return lhs.s < rhs.s;
}


template<typename Simulator, typename InputIt, typename OutputIt>
OutputIt digitize(InputIt first1, InputIt last1, OutputIt d_first, const Configurator& cfg)
{
  std::vector<GeantHit> hits(first1, last1);
  std::vector<DigiChannel> digi_data;

  static Simulator simu(cfg);
  simu.Simulate(hits, digi_data);

  return std::copy(std::begin(digi_data), std::end(digi_data), d_first);
}

}

#endif
