#include <iostream>

#include "TChain.h"
#include "Ttypes.h"

#include "StarMagField/StarMagField.h"
#include "StEvent/StTpcRawData.h"
#include "StTpcDb/StTpcDb.h"
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

  StarMagField starMagField;
  StTpcDb::instance()->SetDriftVelocity();
  StTpcDb::instance()->SetTpcRotations();

  TChain trsTreeChain("t", "tpcrs test TTree");
  trsTreeChain.AddFile("geant_event.root");

  GeantEvent* geantEvent_inp = new GeantEvent();
  GeantEvent  geantEvent_out;

  trsTreeChain.SetBranchAddress("b", &geantEvent_inp);

  std::ofstream logFile_inp("test_tpcrs_inp.log");
  std::ofstream logFile_out("test_tpcrs_out.log");

  for (int iRecord = 1; iRecord <= trsTreeChain.GetEntries(); iRecord++)
  {
    trsTreeChain.GetEntry(iRecord - 1);
    geantEvent_inp->Print(logFile_inp);

    geantEvent_out.hits     = geantEvent_inp->hits;
    geantEvent_out.tracks   = geantEvent_inp->tracks;
    geantEvent_out.vertices = geantEvent_inp->vertices;

    std::cout << geantEvent_inp->hits.size() << std::endl;

    St_g2t_tpc_hit* g2t_tpc_hit = new St_g2t_tpc_hit("g2t_tpc_hit");
    for (auto hit : geantEvent_inp->hits)
      g2t_tpc_hit->AddAt(&hit);

    St_g2t_track* g2t_track = new St_g2t_track("g2t_track");
    for (auto track : geantEvent_inp->tracks)
      g2t_track->AddAt(&track);

    St_g2t_vertex* g2t_vertex = new St_g2t_vertex("g2t_vertex");
    for (auto vertex : geantEvent_inp->vertices)
      g2t_vertex->AddAt(&vertex);

    tpcrs->InitRun(0);
    tpcrs->Make(g2t_tpc_hit, g2t_track, g2t_vertex, tpcraw);

    geantEvent_out.Fill(tpcraw);
    geantEvent_out.Print(logFile_out);
  }

  delete geantEvent_inp;

  return EXIT_SUCCESS;
}
