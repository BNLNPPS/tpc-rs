#include <iostream>

#include "TChain.h"
#include "Ttypes.h"

#include "StEvent/StTpcRawData.h"
#include "StTpcRSMaker/StTpcRSMaker.h"
#include "tables/St_g2t_tpc_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_trigDetSums_Table.h"
#include "tables/St_TpcAvgPowerSupply_Table.h"
#include "tables/St_tpcCorrection_Table.h"
#include "tests/GeantEvent.h"


int main(int argc, char **argv)
{
  StTpcRSMaker* tpcrs = new StTpcRSMaker();
  StTpcRawData* tpcraw = new StTpcRawData();

  TChain trsTreeChain("t", "tpcrs test TTree");
  trsTreeChain.AddFile("geant_event.root");

  GeantEvent* geantEvent = new GeantEvent();

  trsTreeChain.SetBranchAddress("b", &geantEvent);

  std::ofstream logFile("test_tpcrs.log");

  for (int iRecord = 1; iRecord <= trsTreeChain.GetEntries(); iRecord++)
  {
    trsTreeChain.GetEntry(iRecord - 1);
    geantEvent->Print(logFile);

    std::cout << geantEvent->hits.size() << std::endl;

    St_g2t_tpc_hit* g2t_tpc_hit = new St_g2t_tpc_hit("g2t_tpc_hit");
    for (auto hit : geantEvent->hits)
      g2t_tpc_hit->AddAt(&hit);

    St_g2t_track* g2t_track = new St_g2t_track("g2t_track");
    for (auto track : geantEvent->tracks)
      g2t_track->AddAt(&track);

    St_g2t_vertex* g2t_vertex = new St_g2t_vertex("g2t_vertex");
    for (auto vertex : geantEvent->vertices)
      g2t_vertex->AddAt(&vertex);

  }

  delete geantEvent;

  return EXIT_SUCCESS;
}
