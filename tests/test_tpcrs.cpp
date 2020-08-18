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

  std::cout << "test_name:   " << test_name << "\n"
            << "max_records: " << max_records << "\n";

  tpcrs::Configurator cfg(test_name);
  tpcrs::Simulator simulator(cfg);
  tpcrs::MagField mag_field(cfg);

  TChain trsTreeChain("t", "tpcrs test TTree");
  trsTreeChain.AddFile( cfg.Locate(test_name + ".root").c_str() );

  GeantEvent* geantEvent_inp = new GeantEvent();
  GeantEvent  geantEvent_out;

  trsTreeChain.SetBranchAddress("b", &geantEvent_inp);

  std::ofstream logFile_inp(test_name + "_inp.log");
  std::ofstream logFile_out(test_name + "_out.log");

  max_records = (max_records < 0 || max_records > trsTreeChain.GetEntries() ? trsTreeChain.GetEntries() : max_records);

  for (int iRecord = 1; iRecord <= max_records; iRecord++)
  {
    logFile_inp << "event: " << iRecord << "\n";
    logFile_out << "event: " << iRecord << "\n";

    trsTreeChain.GetEntry(iRecord - 1);
    geantEvent_inp->Print(logFile_inp);

    geantEvent_out.hits     = geantEvent_inp->hits;
    geantEvent_out.tracks   = geantEvent_inp->tracks;
    geantEvent_out.vertices = geantEvent_inp->vertices;

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

    geantEvent_out.Fill(digi_data);
    geantEvent_out.Print(logFile_out);
  }

  delete geantEvent_inp;

  return EXIT_SUCCESS;
}


tpcrs::SimulatedHit merge(const g2t_tpc_hit& hit, const g2t_track& particle, const g2t_vertex& vertex)
{
  return tpcrs::SimulatedHit{
    hit.track_p,
    particle.ge_pid,
    hit.volume_id,
    {hit.x[0], hit.x[1], hit.x[2]},
    {hit.p[0], hit.p[1], hit.p[2]},
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
