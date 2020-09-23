#ifndef GeantEvent_h
#define GeantEvent_h

#include <vector>
#include "Rtypes.h"
#include "tpcrs/config_type.h"


struct DigitizedHit
{
  unsigned int sector;
  unsigned int row;
  unsigned int pad;
  std::vector<short> ADCs;
  std::vector<short> IDTs;

  DigitizedHit() : sector(0), row(0), pad(0), ADCs(), IDTs() { }

  DigitizedHit(int s, int r, int p, std::vector<short> adcs, std::vector<short> idts) :
    sector(s), row(r), pad(p), ADCs(adcs), IDTs(idts) { }

  ClassDef(DigitizedHit, 2)
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

      for (unsigned int tb = 0; tb < digiHit.ADCs.size(); tb++) {
        if (!digiHit.ADCs[tb] || !digiHit.IDTs[tb]) continue;
        diff.total++;

        T  digi_inp{digiHit.sector, digiHit.row, digiHit.pad, tb+1, digiHit.ADCs[tb], digiHit.IDTs[tb]};
        const T& digi_out = *dd_iter;

        if (digi_inp.channel < digi_out.channel) {
          diff.unmatched++;
          continue;
        }
        else if (digi_out.channel < digi_inp.channel) {
          diff.unmatched++;
          ++dd_iter;
          tb--; // i.e. do not increment before next iteration
          continue;
        }
        else {
          if (digi_inp.adc != digi_out.adc || digi_inp.track_id != digi_out.track_id) {
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
