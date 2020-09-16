#pragma once

#include <iomanip>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

#include "tpcrs/configurator.h"

namespace tpcrs {

struct TPC
{
  enum Half { first, second };
};


struct DigiChannel
{
  unsigned int sector : 5, row : 7, pad : 10, timebin : 10;
};


inline bool operator< (const DigiChannel& a, const DigiChannel& b)
{
  return std::tie(a.sector, a.row, a.pad, a.timebin) <
         std::tie(b.sector, b.row, b.pad, b.timebin);
}


inline bool operator== (const DigiChannel& a, const DigiChannel& b)
{
  return std::tie(a.sector, a.row, a.pad, a.timebin) ==
         std::tie(b.sector, b.row, b.pad, b.timebin);
}


/**
 * Associates a packed hardware index with the data.
 *
 *                sector   row    pad   timebin
 * projected max      24    72    182       512
 * bits                5 +   7 +   10 +      10 = 32 = 4 bytes = unsigned int
 * max value          32   128   1024      1024
 */
struct DigiHit
{
  DigiChannel channel;
  short adc;
  short idt;
};


struct DistortedHit
{
  /// Distorted position of the simulated hit at the center of a pad row layer
  /// in the global coordinate system
  double x, y, z;

  /// Distorted particle momentum at point (x,y,z)
  double px, py, pz;

  /// Energy deposited by the particle within pad row boundaries
  double de;

  /// Path length within pad row boundaries
  double ds;

  /// Number of primary electrons
  int np;
};


struct SimulatedCharge
{
  float charge;
  short track_id;
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
  double x, y, z;

  /// Local particle momentum at point x
  double px, py, pz;

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
};


inline std::istream& operator>>(std::istream& is, SimulatedHit& hit)
{
  int p = 6; // precision

  is >> std::fixed >> std::setfill('0')
     >> std::setw(4)      >> hit.volume_id
     >> std::setfill(' ')
     >> std::setw(6)      >> hit.track_id
     >> std::setw(p + 6)  >> hit.s
     >> std::setw(6)      >> hit.particle_id
     >> std::setw(p + 6)  >> hit.x
     >> std::setw(p + 6)  >> hit.y
     >> std::setw(p + 6)  >> hit.z
     >> std::setw(p + 6)  >> hit.px
     >> std::setw(p + 6)  >> hit.py
     >> std::setw(p + 6)  >> hit.pz
     >> std::scientific
     >> std::setw(p + 10) >> hit.de
     >> std::fixed
     >> std::setw(p + 6)  >> hit.ds
     >> std::scientific
     >> std::setw(p + 10) >> hit.tof
     >> std::fixed
     >> std::setw(p + 6)  >> hit.lgam;

  return is;
}


inline std::ostream& operator<<(std::ostream& os, const SimulatedHit& hit)
{
  int p = 6; // precision

  os << std::fixed << std::setfill('0')
     << std::setw(4)      << hit.volume_id%10000
     << std::setfill(' ')
     << std::setw(6)      << hit.track_id
     << std::setw(p + 6)  << hit.s
     << std::setw(6)      << hit.particle_id
     << std::setw(p + 6)  << hit.x
     << std::setw(p + 6)  << hit.y
     << std::setw(p + 6)  << hit.z
     << std::setw(p + 6)  << hit.px
     << std::setw(p + 6)  << hit.py
     << std::setw(p + 6)  << hit.pz
     << std::scientific
     << std::setw(p + 10) << hit.de
     << std::fixed
     << std::setw(p + 6)  << hit.ds
     << std::scientific
     << std::setw(p + 10) << hit.tof
     << std::fixed
     << std::setw(p + 6)  << hit.lgam;

  return os;
}


/**
 * Simulated hits must be sorted by sector, track id, and track length up to
 * this hit
 */
inline bool operator< (const SimulatedHit& lhs, const SimulatedHit& rhs)
{
  int lhs_sector = (lhs.volume_id % 10000) / 100;
  int rhs_sector = (rhs.volume_id % 10000) / 100;
  return std::tie(lhs_sector, lhs.track_id, lhs.s) < std::tie(rhs_sector, rhs.track_id, rhs.s);
}


inline SimulatedCharge& operator+= (SimulatedCharge& a, const SimulatedCharge& b)
{
  a.charge += b.charge;
  if (a.track_id != b.track_id && a.charge < 2 * b.charge)
    a.track_id = b.track_id;
  return a;
}


template<typename InputIt, typename OutputIt, typename Simulator>
OutputIt distort(InputIt first_hit, InputIt last_hit, OutputIt d_first, const Simulator& simu)
{
  return simu.Distort(first_hit, last_hit, d_first);
}


template<typename InputIt, typename OutputIt, typename Simulator>
OutputIt digitize(InputIt first_hit, InputIt last_hit, OutputIt d_first, const Simulator& simu)
{
  return simu.Digitize(first_hit, last_hit, d_first);
}


struct DigiChannelMap
{
  DigiChannelMap(const Configurator& cfg, DigiChannel first = {1, 1, 1, 0}, DigiChannel last = {2, 1, 1, 0}) :
    channels(),
    n_sectors(cfg.S<tpcDimensions>().numberOfSectors),
    n_rows(cfg.S<tpcPadPlanes>().innerPadRows + cfg.S<tpcPadPlanes>().outerPadRows),
    n_timebins(cfg.S<tpcElectronics>().numberOfTimeBins),
    num_pads_(),
    num_pads_cumul_(n_rows + 1, 0)
  {
    num_pads_.insert(std::end(num_pads_), cfg.S<tpcPadPlanes>().innerPadsPerRow,
                                          cfg.S<tpcPadPlanes>().innerPadsPerRow + cfg.S<tpcPadPlanes>().innerPadRows),
    num_pads_.insert(std::end(num_pads_), cfg.S<tpcPadPlanes>().outerPadsPerRow,
                                          cfg.S<tpcPadPlanes>().outerPadsPerRow + cfg.S<tpcPadPlanes>().outerPadRows);

    std::partial_sum(std::begin(num_pads_), std::end(num_pads_), std::begin(num_pads_cumul_)+1);

    DigiChannel channel = first;
    while ( !(channel == last) )
    {
      channels.push_back(channel);

      channel.timebin++;

      if (channel.timebin > n_timebins - 1) { channel.timebin = 0; channel.pad++; }
      if (channel.pad > num_pads_[channel.row-1]) { channel.pad = 1; channel.row++; }
      if (channel.row > n_rows) { channel.row = 1; channel.sector++; }
    }
  }

  int n_pads(int row)     const { return num_pads_[row-1]; }
  int total_pads(int row) const { return num_pads_cumul_[row-1]; }
  int total_timebins()    const { return num_pads_cumul_[n_rows] * n_timebins; }

  int n_sectors;
  int n_rows;
  int n_timebins;

  std::vector<DigiChannel> channels;

 private:

  std::vector<int> num_pads_;
  std::vector<int> num_pads_cumul_;
};

}
