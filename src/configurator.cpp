#include <iostream>
#include <sys/stat.h>

#include "tpcrs/configurator.h"
#include "config.h"

namespace tpcrs {

Configurator::Configurator(std::string configname, std::string filename) :
  name(configname)
{
  std::string paths(TPCRS_CONFIG_SEARCH_PATHS);

  size_t last = 0;
  size_t next = 0;
  while ((next = paths.find(':', last)) != std::string::npos)
  {
    search_paths.push_back(paths.substr(last, next-last));
    last = next + 1;
  }
  search_paths.push_back(paths.substr(last));

  if (filename.empty())
    filename = Locate(name + ".yaml");

  yaml = YAML::LoadFile(filename);
}


std::string Configurator::Locate(std::string filename) const
{
  struct stat buffer;   
  for (std::string path : search_paths)
  {
    std::string fullname(path + "/" + filename);
    if (stat(fullname.c_str(), &buffer) == 0)
      return fullname;
  }

  return "";
}


template<> std::string ConfigNodeName<trigDetSums>()          { return "Calibrations/rich/trigDetSums"; };
template<> std::string ConfigNodeName<asic_thresholds>()      { return "Calibrations/tpc/asic_thresholds"; };
template<> std::string ConfigNodeName<ResponseSimulator>()    { return "Calibrations/tpc/ResponseSimulator"; };
template<> std::string ConfigNodeName<tpcAltroParams>()       { return "Calibrations/tpc/tpcAltroParams"; };
template<> std::string ConfigNodeName<tpcAnodeHVavg>()        { return "Calibrations/tpc/tpcAnodeHVavg"; };
template<> std::string ConfigNodeName<tpcCalibResolutions>()  { return "Calibrations/tpc/tpcCalibResolutions"; };
template<> std::string ConfigNodeName<tpcDriftVelocity>()     { return "Calibrations/tpc/tpcDriftVelocity"; };
template<> std::string ConfigNodeName<tpcEffectiveGeom>()     { return "Calibrations/tpc/tpcEffectiveGeom"; };
template<> std::string ConfigNodeName<tpcElectronics>()       { return "Calibrations/tpc/tpcElectronics"; };
template<> std::string ConfigNodeName<tpcGas>()               { return "Calibrations/tpc/tpcGas"; };
template<> std::string ConfigNodeName<tpcHighVoltages>()      { return "Calibrations/tpc/tpcHighVoltages"; };
template<> std::string ConfigNodeName<tpcOmegaTau>()          { return "Calibrations/tpc/tpcOmegaTau"; };
template<> std::string ConfigNodeName<tpcPadGainT0>()         { return "Calibrations/tpc/tpcPadGainT0"; };
template<> std::string ConfigNodeName<tpcPadrowT0>()          { return "Calibrations/tpc/tpcPadrowT0"; };
template<> std::string ConfigNodeName<tpcSectorT0offset>()    { return "Calibrations/tpc/tpcSectorT0offset"; };
template<> std::string ConfigNodeName<TpcResponseSimulator>() { return "Calibrations/tpc/TpcResponseSimulator"; };
template<> std::string ConfigNodeName<trgTimeOffset>()        { return "Conditions/trg/trgTimeOffset"; };
template<> std::string ConfigNodeName<tpcDimensions>()        { return "Geometry/tpc/tpcDimensions"; };
template<> std::string ConfigNodeName<tpcGlobalPosition>()    { return "Geometry/tpc/tpcGlobalPosition"; };
template<> std::string ConfigNodeName<tpcPadPlanes>()         { return "Geometry/tpc/tpcPadPlanes"; };
template<> std::string ConfigNodeName<tpcWirePlanes>()        { return "Geometry/tpc/tpcWirePlanes"; };
template<> std::string ConfigNodeName<MagFactor>()            { return "RunLog/MagFactor"; };
template<> std::string ConfigNodeName<starClockOnl>()         { return "RunLog/onl/starClockOnl"; };
template<> std::string ConfigNodeName<tss_tsspar>()           { return "tpc/tsspars/tsspar"; };


}
