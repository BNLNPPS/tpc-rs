#include <iostream>

#include "TChain.h"
#include "Ttypes.h"
#include "TBenchmark.h"

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
#include "MagFactor.h"

#include "tests/GeantEvent.h"

#include "tpcrs/configurator.h"
#include "tpcrs/config_yaml.h"


int main(int argc, char **argv)
{
  // Process 1st optional argument
  std::string testName(argc > 1 ? argv[1] : "rcf16000_1_100evts.fzd");

  // Process 2nd optional argument
  int maxRecords(argc > 2 ? std::atoi(argv[2]) : -1);

  // Process 3rd optional argument
  double eCutOff(argc > 3 ? std::atof(argv[3]) : 1e-3);

  std::cout << "testName:   " << testName << "\n"
            << "maxRecords: " << maxRecords << "\n"
            << "eCutOff:    " << eCutOff << "\n";

  tpcrs::Configurator::Instance().Configure(testName);

  gBenchmark = new TBenchmark();

  MagFactor_st mf = tpcrs::Configurator::YAML("RunLog/MagFactor").as<MagFactor_st>();
  StarMagField starMagField(StarMagField::kMapped, mf.ScaleFactor);
  StTpcDb::instance()->SetDriftVelocity();
  StTpcDb::instance()->SetTpcRotations();

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

    St_g2t_tpc_hit* g2t_tpc_hit = new St_g2t_tpc_hit("g2t_tpc_hit");
    for (auto hit : geantEvent_inp->hits)
      g2t_tpc_hit->AddAt(&hit);

    St_g2t_track* g2t_track = new St_g2t_track("g2t_track");
    for (auto track : geantEvent_inp->tracks)
      g2t_track->AddAt(&track);

    St_g2t_vertex* g2t_vertex = new St_g2t_vertex("g2t_vertex");
    for (auto vertex : geantEvent_inp->vertices)
      g2t_vertex->AddAt(&vertex);

    StTpcRawData tpcraw;

    tpcrs.Make(g2t_tpc_hit, g2t_track, g2t_vertex, &tpcraw);

    geantEvent_out.Fill(tpcraw);
    geantEvent_out.Print(logFile_out);
  }

  delete geantEvent_inp;

  return EXIT_SUCCESS;
}
