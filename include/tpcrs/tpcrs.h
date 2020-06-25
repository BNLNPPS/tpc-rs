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
  int np;         /* no. of primary electrons */
  double de;      /* energy deposited at hit */
  double ds;      /* path length within pad row */
  float adc;
  float pad;
  float timebin;
  //DigiChannel ch;
};


struct GeantHit
{
  int track_id;   /* Parent track */
  int particle_id;    /* GEANT particle id */
  int volume_id; /* Volume id packed as SSRR, SS = sector, RR = pad row */
  float x[3];    /* Hit center */
  float p[3];    /* Local momentum */
  float de;      /* energy deposited at hit */
  float ds;      /* path length within pad row */
  float len;     /* track length up to this hit */
  double tof;     /// Time of flight including the GEANT vertex production time
  float lgam;     /// Deprecated. ALOG10(GEKin/AMass)
  SimuHit digi;   /// Deprecated.
};


struct TPC
{
  enum Half { first, second };
};


template<typename Simulator, typename InputIt, typename OutputIt>
OutputIt simulate(InputIt first1, InputIt last1, OutputIt d_first, const Configurator& cfg)
{
  static Simulator simu(cfg);

  std::vector<GeantHit> hits(first1, last1);
  std::vector<DigiChannel> digi_data;

  simu.Simulate(hits, digi_data);

  return std::copy(std::begin(digi_data), std::end(digi_data), d_first);
}

}

#endif
