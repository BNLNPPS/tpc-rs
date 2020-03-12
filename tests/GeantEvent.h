#ifndef GeantEvent_h
#define GeantEvent_h

#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

#include "Rtypes.h"

#include "StEvent/StTpcRawData.h"

#include "g2t_tpc_hit.h"
#include "g2t_track.h"
#include "g2t_vertex.h"

#define __MaxNumberOfTimeBins__ 512


struct DigitizedHit
{
  int sector;
  int row;
  int pad;
  std::vector<Short_t>  ADCs;
  std::vector<UShort_t> IDTs;

  DigitizedHit() : sector(0), row(0), pad(0), ADCs(), IDTs() { }

  DigitizedHit(int s, int r, int p, std::vector<Short_t> adcs, std::vector<UShort_t> idts) :
    sector(s), row(r), pad(p), ADCs(adcs), IDTs(idts) { }

  void Print(std::ostream &os, bool resetIndex=false)
  {
    static int index = 0;

    index = resetIndex ? 0 : index;

    // Print only when hit has values
    std::ostringstream oss;
    for (int i = 0; i < ADCs.size(); i++) {
      if (!ADCs[i] || !IDTs[i]) continue;
      oss << "\t" << i << "/" << ADCs[i] << "/" << IDTs[i];
    }

    if (!oss.str().empty())
      os << index  << "\ts/r/p: " << sector << "/" << row << "/" << pad
         << ":" << oss.str() << "\n";

    index++;
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
    Clear();

    for (int i = 0; i < nHits; i++, hit++) {
      hits.push_back(*hit);
    }

    for (int i = 0; i < nTracks; i++, track++) {
      tracks.push_back(*track);
    }

    for (int i = 0; i < nVertices; i++, vertex++) {
      vertices.push_back(*vertex);
    }
  }

  void Fill(StTpcRawData& tpcrd)
  {
    for (int sector = 1; sector <= tpcrd.getNoSectors(); sector++)
    {
      StTpcDigitalSector *digitalSector = tpcrd.GetSector(sector);

      if (!digitalSector) continue;

      int nRows = St_tpcPadConfigC::instance()->numberOfRows(sector);
      for (int row = 1; row <= nRows; row++)
      {
        Int_t nPads = digitalSector->numberOfPadsInRow(row);
        for(Int_t pad = 1; pad <= nPads; pad++)
        {
          static  Short_t ADCs[__MaxNumberOfTimeBins__];
          static UShort_t IDTs[__MaxNumberOfTimeBins__];
          digitalSector->getTimeAdc(row, pad, ADCs, IDTs);

          digiHits.push_back( DigitizedHit(sector, row, pad,
            std::vector<Short_t>(ADCs, ADCs + __MaxNumberOfTimeBins__),
            std::vector<UShort_t>(IDTs, IDTs + __MaxNumberOfTimeBins__) ) );
        }
      }
    }
  }

  void Print(std::ostream &os)
  {
    static int index = 0;

    os << "event: " << index++ << "\n";
    os << std::setprecision(std::numeric_limits<long double>::digits10 + 2);

    for (auto hit : hits) {
      os << hit.id << ", " << hit.x[0] << ", " << hit.x[1] << ", " << hit.x[2] << ", " << hit.de << ", " << hit.ds << "\n";
    }

    bool resetIndex = true;
    for (auto digiHit : digiHits) {
      digiHit.Print(os, resetIndex);
      resetIndex = false;
    }

    os << std::setprecision(6);
  }

  ClassDef(GeantEvent, 1)
};

#endif
