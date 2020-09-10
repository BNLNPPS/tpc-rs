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

  ClassDef(GeantEvent, 1)
};

#endif
