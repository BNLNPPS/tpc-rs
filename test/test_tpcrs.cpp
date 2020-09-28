#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "TChain.h"

#include "GeantEvent.h"
#include "tpcrs/tpcrs.h"


tpcrs::SimulatedHit merge(const g2t_tpc_hit& hit, const g2t_track& particle, const g2t_vertex& vertex);
tpcrs::SimulatedHit merge(const g2t_tpc_hit& hit, const std::vector<g2t_track>& particles, const std::vector<g2t_vertex>& vertices);


int main(int argc, char **argv)
{
  // Process 1st optional argument
  std::string test_name(argc > 1 ? argv[1] : "starY16_dAu200");

  // Process 2nd optional argument
  int max_records(argc > 2 ? std::atoi(argv[2]) : -1);

  tpcrs::Configurator cfg(test_name);
  tpcrs::Simulator simulator(cfg);
  tpcrs::MagField mag_field(cfg);

  TChain trsTreeChain("t", "tpcrs test TTree");
  trsTreeChain.AddFile( cfg.Locate(test_name + ".root").c_str() );

  GeantEvent* geantEvent_inp = new GeantEvent();

  trsTreeChain.SetBranchAddress("b", &geantEvent_inp);

  max_records = (max_records < 0 || max_records > trsTreeChain.GetEntries() ? trsTreeChain.GetEntries() : max_records);

  std::cout << "test name: " << test_name << "\n"
            << "records: " << max_records << " of " << trsTreeChain.GetEntries() << "\n";

  struct Diff {
    int total;
    int unmatched;
    Diff& operator +=(const Diff& b) { total += b.total; unmatched += b.unmatched; }
  } counts{};

  for (int iRecord = 1; iRecord <= max_records; iRecord++)
  {
    trsTreeChain.GetEntry(iRecord - 1);

    auto convert = [geantEvent_inp](const g2t_tpc_hit& hit) -> tpcrs::SimulatedHit
    {
      return merge(hit, geantEvent_inp->tracks, geantEvent_inp->vertices);
    };

    std::vector<tpcrs::SimulatedHit> hits;
    std::transform(begin(geantEvent_inp->hits), end(geantEvent_inp->hits), std::back_inserter(hits), convert);

    // Sort generated hits
    std::stable_sort(begin(hits), end(hits));

    std::vector<tpcrs::DigiHit>  digi_data;
    simulator.Digitize(std::begin(hits), std::end(hits), back_inserter(digi_data), mag_field);

    counts += geantEvent_inp->Diff<Diff>(digi_data);
  }

  delete geantEvent_inp;

  std::cout << "hits total: " << counts.total << '\n'
            << "hits unmatched: " << counts.unmatched << '\n';

  return counts.unmatched;
}


tpcrs::SimulatedHit merge(const g2t_tpc_hit& hit, const g2t_track& particle, const g2t_vertex& vertex)
{
  return tpcrs::SimulatedHit{
    hit.track_p,
    particle.ge_pid,
    hit.volume_id,
    hit.x[0], hit.x[1], hit.x[2],
    hit.p[0], hit.p[1], hit.p[2],
    hit.de,
    hit.ds,
    hit.length,
    double(hit.tof) + double(vertex.ge_tof),
    hit.lgam
  };
}


tpcrs::SimulatedHit merge(const g2t_tpc_hit& hit, const std::vector<g2t_track>& particles, const std::vector<g2t_vertex>& vertices)
{
  int particle_idx  = hit.track_p;
  int vertex_idx = particles[particle_idx - 1].start_vertex_p;

  return merge(hit, particles[particle_idx - 1], vertices[vertex_idx - 1]);
}
