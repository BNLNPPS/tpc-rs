#ifndef GeantEvent_h
#define GeantEvent_h

#include <cstring>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <vector>

#include "Rtypes.h"

#include "tpcrs/digi_data.h"
#include "tpcrs/structs.h"

#define __MaxNumberOfTimeBins__ 512


struct DigitizedHit
{
  int sector;
  int row;
  int pad;
  std::vector<short>  ADCs;
  std::vector<unsigned short> IDTs;

  DigitizedHit() : sector(0), row(0), pad(0), ADCs(), IDTs() { }

  DigitizedHit(int s, int r, int p, std::vector<short> adcs, std::vector<unsigned short> idts) :
    sector(s), row(r), pad(p), ADCs(adcs), IDTs(idts) { }

  void Print(std::ostream &os)
  {
    // Print only when hit has values
    std::ostringstream oss;
    for (int i = 0; i < ADCs.size(); i++) {
      if (!ADCs[i] || !IDTs[i]) continue;
      oss << "\t" << i << "/" << ADCs[i] << "/" << IDTs[i];
    }

    if (!oss.str().empty())
      os << "s/r/p: " << sector << "/" << row << "/" << pad
         << ":" << oss.str() << "\n";
  }

  ClassDef(DigitizedHit, 1)
};


struct GeantEvent
{
  std::vector<g2t_tpc_hit_st> hits;
  std::vector<g2t_track_st> tracks;
  std::vector<g2t_vertex_st> vertices;

  std::vector<DigitizedHit> digiHits;

  void Clear() {
    hits.clear();
    tracks.clear();
    vertices.clear();
    digiHits.clear();
  }

  void Fill(const g2t_tpc_hit_st* hit, int nHits, const g2t_track_st* track, int nTracks, const g2t_vertex_st* vertex, int nVertices)
  {
    hits.clear();
    for (int i = 0; i < nHits; i++, hit++)
      hits.push_back(*hit);

    tracks.clear();
    for (int i = 0; i < nTracks; i++, track++)
      tracks.push_back(*track);

    vertices.clear();
    for (int i = 0; i < nVertices; i++, vertex++)
      vertices.push_back(*vertex);
  }

  void Fill(tpcrs::DigiData& digi_data)
  {
    digiHits.clear();

    static  short ADCs[__MaxNumberOfTimeBins__];
    static unsigned short IDTs[__MaxNumberOfTimeBins__];

    for (auto tc = digi_data.channels().begin(); tc != digi_data.channels().end(); ++tc)
    {
      ADCs[tc->timebin] = tc->adc;
      IDTs[tc->timebin] = tc->idt;

      auto next_tc = std::next(tc);
      bool same_pad = (tc->sector == next_tc->sector && tc->row == next_tc->row && tc->pad == next_tc->pad);

      if (!same_pad) {
        digiHits.push_back( DigitizedHit(tc->sector, tc->row, tc->pad,
            std::vector<short>(ADCs, ADCs + __MaxNumberOfTimeBins__),
            std::vector<unsigned short>(IDTs, IDTs + __MaxNumberOfTimeBins__) ) );
        std::memset(ADCs, 0, __MaxNumberOfTimeBins__*sizeof(short));
        std::memset(IDTs, 0, __MaxNumberOfTimeBins__*sizeof(unsigned short));
      }
    }
  }

  void Print(std::ostream &os)
  {
    os << std::setprecision(std::numeric_limits<long double>::digits10 + 2);

    for (auto hit : hits) {
      os << hit.id << ", " << hit.x[0] << ", " << hit.x[1] << ", " << hit.x[2] << ", " << hit.de << ", " << hit.ds << "\n";
    }

    for (auto digiHit : digiHits) {
      digiHit.Print(os);
    }

    os << std::setprecision(6);
  }

  ClassDef(GeantEvent, 1)
};

#endif
