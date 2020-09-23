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

    auto dd_iter = begin(digi_data);

    for (auto digiHit : digiHits) {

      for (int i = 0; i < digiHit.ADCs.size(); i++) {
        if (!digiHit.ADCs[i] || !digiHit.IDTs[i]) continue;
        diff.total++;

        T  digi_inp{digiHit.sector, digiHit.row, digiHit.pad, i, digiHit.ADCs[i], digiHit.IDTs[i]};
        const T& digi_out = *dd_iter;

        if (digi_inp.channel < digi_out.channel) {
          diff.unmatched++;
          continue;
        }
        else if (digi_out.channel < digi_inp.channel) {
          diff.unmatched++;
          ++dd_iter;
          i--;
          continue;
        }
        else {
          if (digi_inp.adc != digi_out.adc || digi_inp.idt != digi_out.idt) {
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
