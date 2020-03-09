#include <iostream>

#include "TInterpreter.h"
#include "TSystem.h"

#include "StChain/StMaker.h"
#include "St_db_Maker/StValiSet.h"


StChain* StMaker::gStChain = new StChain();


TTable* StChain::GetDataBase(std::string path, const TDatime *td)
{
  std::string paths("./StarDb/:../StarDb/:/home/smirnovd/tpc-rs/StarDb/");

  std::string sfile( gSystem->Which(paths.c_str(), (path+".C").c_str(), kReadPermission) );

  std::string command(".L " + sfile);

  gInterpreter->ProcessLine(command.c_str());

  TDataSet* ds = reinterpret_cast<TDataSet*>(gInterpreter->Calc("CreateTable()"));

  StValiSet *vs = new StValiSet(("." + std::string(gSystem->BaseName(path.c_str()))).c_str(), ds);
  vs->fTimeMin = TDatime(kMinTime, 0);
  vs->fTimeMax = TDatime(kMaxTime, 0);
  vs->fVers++;

  ds->Add(vs);

  return reinterpret_cast<TTable*>(ds);
}
