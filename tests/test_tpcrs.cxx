#include <iostream>

#include "TChain.h"

#include "StTpcRSMaker/StTpcRSMaker.h"
#include "tests/GeantEvent.h"


int main(int argc, char **argv)
{
  StTpcRSMaker* tpcrs = new StTpcRSMaker();

  TChain trsTreeChain("t", "tpcrs test TTree");
  trsTreeChain.AddFile("geant_event.root");

  GeantEvent* geantEvent = new GeantEvent();

  trsTreeChain.SetBranchAddress("b", &geantEvent);

  const int nRecords = trsTreeChain.GetEntries();

  std::cout << "nRecords: " << nRecords << std::endl;

  std::ofstream logFile("test_tpcrs.log");

  for (int iRecord = 1; iRecord <= nRecords; iRecord++)
  {
    trsTreeChain.GetEntry(iRecord - 1);
    geantEvent->Print(logFile);
  }

  delete geantEvent;

  return EXIT_SUCCESS;
}
