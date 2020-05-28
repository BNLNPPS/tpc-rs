#ifndef tpcrs_configurator_h
#define tpcrs_configurator_h

#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"
#include "tpcrs/structs.h"


namespace tpcrs {

template<typename Struct> std::string ConfigNodeName() { return "undefined"; };

template<> std::string ConfigNodeName<trigDetSums>();
template<> std::string ConfigNodeName<asic_thresholds>();
template<> std::string ConfigNodeName<ResponseSimulator>();
template<> std::string ConfigNodeName<tpcAltroParams>();
template<> std::string ConfigNodeName<tpcDriftVelocity>();
template<> std::string ConfigNodeName<tpcEffectiveGeom>();
template<> std::string ConfigNodeName<tpcElectronics>();
template<> std::string ConfigNodeName<tpcGas>();
template<> std::string ConfigNodeName<tpcPadrowT0>();
template<> std::string ConfigNodeName<TpcResponseSimulator>();
template<> std::string ConfigNodeName<tpcDimensions>();
template<> std::string ConfigNodeName<tpcPadPlanes>();
template<> std::string ConfigNodeName<tpcWirePlanes>();
template<> std::string ConfigNodeName<MagFactor>();
template<> std::string ConfigNodeName<starClockOnl>();
template<> std::string ConfigNodeName<tss_tsspar>();
template<> std::string ConfigNodeName<iTPCSurvey>();

/**
 * A singleton class 
 */
class Configurator
{
 public:

  Configurator(std::string cfgname = "");

  std::string Locate(std::string filename = "") const;

  std::string File() const { return Locate(name + ".yaml"); }

  YAML::Node YAML(std::string taxon) const { return yaml[taxon]; }

  template<typename Chair>
  Chair& C() const
  {
    static Chair c(*this);

    if (this != &c.cfg_)
      c = Chair(*this);

    return c;
  }

  template<typename Struct>
  const Struct& S(int i = 0) const
  {
    static std::vector<Struct> rows;
    static std::string name{ConfigNodeName<Struct>()};
  
    if (rows.size()) {
      return rows[i];
    }
    else {
      try {
        rows.push_back( YAML(name).as<Struct>() );
      } catch (std::exception& e) {
        rows = YAML(name).as< std::vector<Struct> >();
      }
    }
  
    return rows[i];
  }

 private:

  /// A unique name used in various file names associated with this Configurator.
  std::string name;

  std::vector<std::string> searchPaths;

  YAML::Node yaml;
};

}

#endif
