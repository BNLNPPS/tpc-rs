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

#include "tpcrs/config_type.h"

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


typedef g2t_tpc_hit g2t_tpc_hit_st;
typedef g2t_track   g2t_track_st;
typedef g2t_vertex  g2t_vertex_st;

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

  template<typename T>
  void Fill(const std::vector<T>& digi_data)
  {
    digiHits.clear();

    static  short ADCs[__MaxNumberOfTimeBins__];
    static unsigned short IDTs[__MaxNumberOfTimeBins__];

    for (auto tc = begin(digi_data); tc != end(digi_data); ++tc)
    {
      ADCs[tc->channel.timebin] = tc->adc;
      IDTs[tc->channel.timebin] = tc->idt;

      auto next_tc = std::next(tc);
      bool same_pad = (tc->channel.sector == next_tc->channel.sector && tc->channel.row == next_tc->channel.row && tc->channel.pad == next_tc->channel.pad);

      if (!same_pad) {
        digiHits.push_back( DigitizedHit(tc->channel.sector, tc->channel.row, tc->channel.pad,
            std::vector<short>(ADCs, ADCs + __MaxNumberOfTimeBins__),
            std::vector<unsigned short>(IDTs, IDTs + __MaxNumberOfTimeBins__) ) );
        std::memset(ADCs, 0, __MaxNumberOfTimeBins__*sizeof(short));
        std::memset(IDTs, 0, __MaxNumberOfTimeBins__*sizeof(unsigned short));
      }
    }
  }

  template<typename D, typename T>
  D Diff(const std::vector<T>& digi_data)
  {
    D diff{};

    using DigiChannelType = decltype((T().channel));

    auto dd_iter = begin(digi_data);

    for (auto digiHit : digiHits) {

      for (int i = 0; i < digiHit.ADCs.size(); i++) {
        if (!digiHit.ADCs[i] || !digiHit.IDTs[i]) continue;
        diff.total++;

        DigiChannelType  channel_inp{digiHit.sector, digiHit.row, digiHit.pad, i};
        const DigiChannelType& channel_out = dd_iter->channel;

        if (channel_inp < channel_out) {
          diff.unmatched++;


          continue;
        }
        else if (channel_out < channel_inp) {
          diff.unmatched++;

          ++dd_iter;
          i--;
          continue;
        }
        else {
          if (digiHit.ADCs[i] != dd_iter->adc || digiHit.IDTs[i] != dd_iter->idt) {
            diff.unmatched++;
          }
        }
        ++dd_iter;
      }
    }

    return diff;
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
