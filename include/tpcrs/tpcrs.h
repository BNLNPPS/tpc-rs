#ifndef TPCRS_TPCRS_H_
#define TPCRS_TPCRS_H_

#include <vector>
#include <tuple>

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


struct SimulatedHit
{
  /// Unique id of the particle produced this hit. The simulated signal is
  /// associated with tracks using these ids
  int track_id;

  /// Physical particle id. It is primarily used to get the particle's charge
  /// and mass from a local list of known particles
  int particle_id;

  /// Volume id packed as a four digit decimal number SSRR, SS = sector, RR = pad row
  /// Many configuration properties are defined on per sector basis. The prior
  /// knowledge of hits' pad row allows to split loopers into multiple shorter
  /// tracks
  int volume_id;

  /// Position of the simulated hit at the center of a pad row layer in the
  /// global coordinate system
  double x[3];

  /// Local particle momentum at point x
  double p[3];

  /// Energy deposited by the particle within pad row boundaries. Not strictly
  /// necessary as it is used only in special cases such as "stopped electrons"
  /// and indicated by negative `de` and `ds` < 50 um
  double de;

  /// Path length within pad row boundaries. Used to define beginning and end
  /// points of the track rebuilt around the distorted hit. Also used in
  /// calculation of local gain correction due to de/ds effects
  double ds;

  /// Distance from the origin to the point x. Not strictly necessary as it is
  /// used only for hit ordering
  double s;

  /// Time of flight up to the point x including the vertex production time. The
  /// time of flight is used to correct the hit `z` position for drift velocity
  /// in cm/us obtained from laser data
  double tof;

  /// log_10(E_kin/mass) -- Deprecated
  float lgam;

  /// Deprecated
  SimuHit digi;
};


struct TPC
{
  enum Half { first, second };
};


/**
 * Expects a type similar to GeneratedHit
 */
inline bool operator< (const SimulatedHit& lhs, const SimulatedHit& rhs)
{
  int lhs_sector = (lhs.volume_id % 10000) / 100;
  int rhs_sector = (rhs.volume_id % 10000) / 100;
  return std::tie(lhs_sector, lhs.track_id, lhs.s) < std::tie(rhs_sector, rhs.track_id, rhs.s);
}


inline bool operator< (const SimulatedCharge& lhs, const SimulatedCharge& rhs)
{
  return std::tie(lhs.sector, lhs.row, lhs.pad, lhs.timebin) < std::tie(rhs.sector, rhs.row, rhs.pad, rhs.timebin);
}


template<typename Simulator, typename InputIt, typename OutputIt>
OutputIt digitize(InputIt first1, InputIt last1, OutputIt d_first, const Configurator& cfg)
{
  static Simulator simu(cfg);
  simu.Simulate(first1, last1, d_first);

  return d_first;
}

}

#endif
