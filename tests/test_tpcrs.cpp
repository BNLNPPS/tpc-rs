#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "TChain.h"

#include "tests/GeantEvent.h"
#include "tpcrs/tpcrs.h"
#include "tpcrs/configurator.h"
#include "tpcrs/tpcrs.h"
#include "simulator.h"


tpcrs::GeantHit merge(const g2t_tpc_hit& hit, const g2t_track& particle, const g2t_vertex& vertex);
tpcrs::GeantHit merge(const g2t_tpc_hit& hit, const std::vector<g2t_track>& particles, const std::vector<g2t_vertex>& vertices);


int main(int argc, char **argv)
{
  // Process 1st optional argument
  std::string testName(argc > 1 ? argv[1] : "starY16_dAu200");

  // Process 2nd optional argument
  int maxRecords(argc > 2 ? std::atoi(argv[2]) : -1);

  // Process 3rd optional argument
  double eCutOff(argc > 3 ? std::atof(argv[3]) : 1e-3);

  std::cout << "testName:   " << testName << "\n"
            << "maxRecords: " << maxRecords << "\n"
            << "eCutOff:    " << eCutOff << "\n";

  tpcrs::Configurator::Instance().Configure(testName);

  StTpcRSMaker tpcrs(eCutOff);

  TChain trsTreeChain("t", "tpcrs test TTree");
  trsTreeChain.AddFile( tpcrs::Configurator::Instance().Locate(testName + ".root").c_str() );

  GeantEvent* geantEvent_inp = new GeantEvent();
  GeantEvent  geantEvent_out;

  trsTreeChain.SetBranchAddress("b", &geantEvent_inp);

  std::ofstream logFile_inp(testName + "_inp.log");
  std::ofstream logFile_out(testName + "_out.log");

  maxRecords = (maxRecords < 0 || maxRecords > trsTreeChain.GetEntries() ? trsTreeChain.GetEntries() : maxRecords);

  for (int iRecord = 1; iRecord <= maxRecords; iRecord++)
  {
    logFile_inp << "event: " << iRecord << "\n";
    logFile_out << "event: " << iRecord << "\n";

    trsTreeChain.GetEntry(iRecord - 1);
    geantEvent_inp->Print(logFile_inp);

    geantEvent_out.hits     = geantEvent_inp->hits;
    geantEvent_out.tracks   = geantEvent_inp->tracks;
    geantEvent_out.vertices = geantEvent_inp->vertices;

    std::cout << geantEvent_inp->hits.size() << "\n";

    auto convert = [geantEvent_inp](const g2t_tpc_hit& hit) -> tpcrs::GeantHit
    {
      return merge(hit, geantEvent_inp->tracks, geantEvent_inp->vertices);
    };

    std::vector<tpcrs::GeantHit> hits;
    std::transform(begin(geantEvent_inp->hits), end(geantEvent_inp->hits), std::back_inserter(hits), convert);

    tpcrs::DigiData  digidata;
    tpcrs.Make(hits, digidata);

    geantEvent_out.Fill(digidata);
    geantEvent_out.Print(logFile_out);
  }

  delete geantEvent_inp;

  return EXIT_SUCCESS;
}


tpcrs::GeantHit merge(const g2t_tpc_hit& hit, const g2t_track& particle, const g2t_vertex& vertex)
{
  return tpcrs::GeantHit{
    hit.track_p,
    particle.ge_pid,
    hit.volume_id,
    {hit.x[0], hit.x[1], hit.x[2]},
    {hit.p[0], hit.p[1], hit.p[2]},
    hit.de,
    hit.ds,
    hit.length,
    hit.tof + vertex.ge_tof,
    hit.lgam,
    0,         // OUTPUT: adc
    0,         // OUTPUT: pad
    0,         // OUTPUT: timebin
    0          // OUTPUT: np
  };
}


tpcrs::GeantHit merge(const g2t_tpc_hit& hit, const std::vector<g2t_track>& particles, const std::vector<g2t_vertex>& vertices)
{
  int particle_idx  = hit.track_p;
  int vertex_idx = particles[particle_idx - 1].start_vertex_p;

  return merge(hit, particles[particle_idx - 1], vertices[vertex_idx - 1]);
}
