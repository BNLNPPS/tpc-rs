#include <iostream>


#include "TSystem.h"
#include "TTable.h"

#include "StChain/StMaker.h"
#include "St_db_Maker/StValiSet.h"

#include "config_structs.h"
#include "tpcrs/configurator.h"


StChain* StMaker::gStChain = new StChain();


TTable* StChain::GetDataBase(std::string path, const TDatime *td)
{
  static YAML::Node config = YAML::LoadFile(tpcrs::Configurator::File());

  TTable* table;

  for (auto cst : tpcrs::configStructs)
  {
    if (path != cst->tname) continue;
    table = cst->yaml2table(tpcrs::Configurator::YAML(cst->tname));
    break;
  }

  TDataSet* ds = reinterpret_cast<TDataSet*>(table);

  StValiSet *vs = new StValiSet(("." + std::string(gSystem->BaseName(path.c_str()))).c_str(), ds);
  vs->fTimeMin = TDatime(kMinTime, 0);
  vs->fTimeMax = TDatime(kMaxTime, 0);
  vs->fVers++;

  ds->Add(vs);

  return reinterpret_cast<TTable*>(ds);
}
