#include <iostream>

#include "TInterpreter.h"
#include "TSystem.h"

#include "StChain/StMaker.h"


StChain* StMaker::gStChain = new StChain();


TTable* StChain::GetDataBase(std::string path, const TDatime *td)
{
  std::string paths("./StarDb/:../StarDb/:/home/smirnovd/tpc-rs/StarDb/");

  std::string sfile( gSystem->Which(paths.c_str(), (path+".C").c_str(), kReadPermission) );

  std::string command(".L " + sfile);

  gInterpreter->ProcessLine(command.c_str());

  return reinterpret_cast<TTable*>(gInterpreter->Calc("CreateTable()"));
}
