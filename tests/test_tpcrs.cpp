#include <fstream>
#include <iostream>
#include <string>

#include "TChain.h"

#include "tests/GeantEvent.h"
#include "tpcrs/configurator.h"
#include "tpcrs/tpcrs.h"
#include "simulator.h"


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

    tpcrs::DigiData  digidata;
    tpcrs.Make(geantEvent_inp->hits, geantEvent_inp->tracks, geantEvent_inp->vertices, digidata);

    geantEvent_out.Fill(digidata);
    geantEvent_out.Print(logFile_out);
  }

  delete geantEvent_inp;

  return EXIT_SUCCESS;
}
