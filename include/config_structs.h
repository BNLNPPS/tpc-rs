#ifndef tpcrs_config_structs_h
#define tpcrs_config_structs_h

#include <exception>
#include <string>
#include <vector>

#include "TTable.h"

#include "yaml-cpp/yaml.h"

#include "tables/St_MDFCorrection_Table.h"
#include "tables/St_MagFactor_Table.h"
#include "tables/St_Survey_Table.h"
#include "tables/St_TpcAvgCurrent_Table.h"
#include "tables/St_TpcAvgPowerSupply_Table.h"
#include "tables/St_TpcEffectivedX_Table.h"
#include "tables/St_TpcResponseSimulator_Table.h"
#include "tables/St_TpcSecRowCor_Table.h"
#include "tables/St_asic_thresholds_Table.h"
#include "tables/St_beamInfo_Table.h"
#include "tables/St_g2t_tpc_hit_Table.h"
#include "tables/St_g2t_track_Table.h"
#include "tables/St_g2t_vertex_Table.h"
#include "tables/St_iTPCSurvey_Table.h"
#include "tables/St_itpcPadGainT0_Table.h"
#include "tables/St_richvoltages_Table.h"
#include "tables/St_spaceChargeCor_Table.h"
#include "tables/St_starClockOnl_Table.h"
#include "tables/St_starMagOnl_Table.h"
#include "tables/St_tpcAltroParams_Table.h"
#include "tables/St_tpcAnodeHV_Table.h"
#include "tables/St_tpcAnodeHVavg_Table.h"
#include "tables/St_tpcCalibResolutions_Table.h"
#include "tables/St_tpcChargeEvent_Table.h"
#include "tables/St_tpcCorrection_Table.h"
#include "tables/St_tpcDimensions_Table.h"
#include "tables/St_tpcDriftVelocity_Table.h"
#include "tables/St_tpcEffectiveGeom_Table.h"
#include "tables/St_tpcElectronics_Table.h"
#include "tables/St_tpcFieldCage_Table.h"
#include "tables/St_tpcFieldCageShort_Table.h"
#include "tables/St_tpcGas_Table.h"
#include "tables/St_tpcGlobalPosition_Table.h"
#include "tables/St_tpcGridLeak_Table.h"
#include "tables/St_tpcHVPlanes_Table.h"
#include "tables/St_tpcHighVoltages_Table.h"
#include "tables/St_tpcOmegaTau_Table.h"
#include "tables/St_tpcPadConfig_Table.h"
#include "tables/St_tpcPadGainT0_Table.h"
#include "tables/St_tpcPadPlanes_Table.h"
#include "tables/St_tpcPadResponse_Table.h"
#include "tables/St_tpcPadrowT0_Table.h"
#include "tables/St_tpcPedestal_Table.h"
#include "tables/St_tpcRDOMap_Table.h"
#include "tables/St_tpcRDOMasks_Table.h"
#include "tables/St_tpcRDOT0offset_Table.h"
#include "tables/St_tpcSCGL_Table.h"
#include "tables/St_tpcSectorPosition_Table.h"
#include "tables/St_tpcSectorT0offset_Table.h"
#include "tables/St_tpcSlowControlSim_Table.h"
#include "tables/St_tpcWirePlanes_Table.h"
#include "tables/St_trgTimeOffset_Table.h"
#include "tables/St_trigDetSums_Table.h"
#include "tables/St_tss_tsspar_Table.h"

#include "tpcrs/config_yaml.h"


namespace tpcrs {

struct ConfigStructI
{
  ConfigStructI(std::string name, std::string name2="") : tname(name), tname2(name2) { }

  virtual bool table2yaml(const TTable& rp, YAML::Node& node) = 0;
  virtual TTable* yaml2table(const YAML::Node& node) = 0;

  std::string tname;
  std::string tname2;
};


template<typename Struct_t, typename Table_t>
struct ConfigStruct : public ConfigStructI
{
  ConfigStruct(std::string name, std::string name2="") : ConfigStructI(name, name2) { }

  virtual bool table2yaml(const TTable& rp, YAML::Node& node)
  {
    if (rp.GetSize() == 1)
    {
      node[tname] = *static_cast<const Struct_t*>(rp.GetArray());
      return true;
    }
    else if ( rp.GetSize() > 1 )
    {
      std::vector<Struct_t> rows;

      for (int i = 0; i < rp.GetSize(); ++i)
        rows.push_back( *static_cast<const Struct_t*>( rp[i] ) );

      node[tname] = rows;
      return true;
    }

    return false;
  }


  virtual TTable* yaml2table(const YAML::Node& node)
  {
    TTable* table = nullptr;
    size_t slash_pos =  tname.find_last_of('/');
    std::string tbase = tname.substr( slash_pos+1 );

    if (!node) {
      table = new Table_t(tbase.c_str(), 0);
      table->Mark();
      return table;
    }

    try
    {
      Struct_t row = node.as< Struct_t >();
      table = new Table_t(tbase.c_str(), 1);
      table->AddAt(&row);
    }
    catch (std::exception& e)
    {
      std::vector<Struct_t> rows = node.as< std::vector<Struct_t> >();
      table = new Table_t(tbase.c_str(), rows.size());
      for (Struct_t& row : rows) table->AddAt(&row);
    }

    return table;
  }

};


std::vector<ConfigStructI*> configStructs{
  new ConfigStruct<richvoltages_st,        St_richvoltages>         ("Calibrations/rich/richvoltages"),
  new ConfigStruct<spaceChargeCor_st,      St_spaceChargeCor>       ("Calibrations/rich/spaceChargeCor"),
  new ConfigStruct<spaceChargeCor_st,      St_spaceChargeCor>       ("Calibrations/rich/spaceChargeCorR2"),
  new ConfigStruct<trigDetSums_st,         St_trigDetSums>          ("Calibrations/rich/trigDetSums"),
  new ConfigStruct<asic_thresholds_st,     St_asic_thresholds>      ("Calibrations/tpc/asic_thresholds"),
  new ConfigStruct<itpcPadGainT0_st,       St_itpcPadGainT0>        ("Calibrations/tpc/itpcPadGainT0"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcAdcCorrectionB"),
  new ConfigStruct<MDFCorrection_st,       St_MDFCorrection>        ("Calibrations/tpc/TpcAdcCorrectionMDF"),
  new ConfigStruct<tpcAltroParams_st,      St_tpcAltroParams>       ("Calibrations/tpc/tpcAltroParams"),
  new ConfigStruct<tpcAnodeHV_st,          St_tpcAnodeHV>           ("Calibrations/tpc/tpcAnodeHV"),
  new ConfigStruct<tpcAnodeHVavg_st,       St_tpcAnodeHVavg>        ("Calibrations/tpc/tpcAnodeHVavg"),
  new ConfigStruct<TpcAvgCurrent_st,       St_TpcAvgCurrent>        ("Calibrations/tpc/TpcAvgCurrent"),
  new ConfigStruct<TpcAvgPowerSupply_st,   St_TpcAvgPowerSupply>    ("Calibrations/tpc/TpcAvgPowerSupply"),
  new ConfigStruct<tpcCalibResolutions_st, St_tpcCalibResolutions>  ("Calibrations/tpc/tpcCalibResolutions"),
  new ConfigStruct<tpcChargeEvent_st,      St_tpcChargeEvent>       ("Calibrations/tpc/tpcChargeEvent"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcCurrentCorrection"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcCurrentCorrectionX"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcdCharge"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcdEdxCor"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcDriftDistOxygen"),
  new ConfigStruct<tpcDriftVelocity_st,    St_tpcDriftVelocity>     ("Calibrations/tpc/tpcDriftVelocity"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcdXCorrectionB"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcEdge"),
  new ConfigStruct<TpcEffectivedX_st,      St_TpcEffectivedX>       ("Calibrations/tpc/TpcEffectivedX"),
  new ConfigStruct<tpcEffectiveGeom_st,    St_tpcEffectiveGeom>     ("Calibrations/tpc/tpcEffectiveGeom", "Calibrations/tpc/tpcEffectiveGeomB"),
  new ConfigStruct<tpcElectronics_st,      St_tpcElectronics>       ("Calibrations/tpc/tpcElectronics", "Calibrations/tpc/tpcElectronicsB"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/tpcGainCorrection"),
  new ConfigStruct<tpcGas_st,              St_tpcGas>               ("Calibrations/tpc/tpcGas"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/tpcGasTemperature"),
  new ConfigStruct<tpcGridLeak_st,         St_tpcGridLeak>          ("Calibrations/tpc/tpcGridLeak"),
  new ConfigStruct<tpcHighVoltages_st,     St_tpcHighVoltages>      ("Calibrations/tpc/tpcHighVoltages"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcLengthCorrectionB"),
  new ConfigStruct<MDFCorrection_st,       St_MDFCorrection>        ("Calibrations/tpc/TpcLengthCorrectionMDF"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/tpcMethaneIn"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcMultiplicity"),
  new ConfigStruct<tpcOmegaTau_st,         St_tpcOmegaTau>          ("Calibrations/tpc/tpcOmegaTau"),
  new ConfigStruct<MDFCorrection_st,       St_MDFCorrection>        ("Calibrations/tpc/TpcPadCorrectionMDF"),
  new ConfigStruct<tpcPadGainT0_st,        St_tpcPadGainT0>         ("Calibrations/tpc/tpcPadGainT0"),
  new ConfigStruct<tpcPadResponse_st,      St_tpcPadResponse>       ("Calibrations/tpc/tpcPadResponse"),
  new ConfigStruct<tpcPadrowT0_st,         St_tpcPadrowT0>          ("Calibrations/tpc/tpcPadrowT0", "Calibrations/tpc/tpcPadrowT0B"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcPhiDirection"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/tpcPressureB"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcrCharge"),
  new ConfigStruct<tpcRDOMap_st,           St_tpcRDOMap>            ("Calibrations/tpc/tpcRDOMap"),
  new ConfigStruct<tpcRDOT0offset_st,      St_tpcRDOT0offset>       ("Calibrations/tpc/tpcRDOT0offset"),
  new ConfigStruct<TpcResponseSimulator_st,St_TpcResponseSimulator> ("Calibrations/tpc/TpcResponseSimulator"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcRowQ"),
  new ConfigStruct<tpcSCGL_st,             St_tpcSCGL>              ("Calibrations/tpc/tpcSCGL"),
  new ConfigStruct<TpcSecRowCor_st,        St_TpcSecRowCor>         ("Calibrations/tpc/TpcSecRowB"),
  new ConfigStruct<TpcSecRowCor_st,        St_TpcSecRowCor>         ("Calibrations/tpc/TpcSecRowC"),
  new ConfigStruct<tpcSectorT0offset_st,   St_tpcSectorT0offset>    ("Calibrations/tpc/tpcSectorT0offset"),
  new ConfigStruct<tpcSlowControlSim_st,   St_tpcSlowControlSim>    ("Calibrations/tpc/tpcSlowControlSim"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcSpaceCharge"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcTanL"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/tpcTimeDependence"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/tpcWaterOut"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcZCorrectionB"),
  new ConfigStruct<tpcCorrection_st,       St_tpcCorrection>        ("Calibrations/tpc/TpcZDC"),
  new ConfigStruct<trgTimeOffset_st,       St_trgTimeOffset>        ("Conditions/trg/trgTimeOffset", "Conditions/trg/trgTimeOffsetB"),
  new ConfigStruct<iTPCSurvey_st,          St_iTPCSurvey>           ("Geometry/tpc/iTPCSurvey"),
  new ConfigStruct<tpcDimensions_st,       St_tpcDimensions>        ("Geometry/tpc/tpcDimensions"),
  new ConfigStruct<tpcFieldCage_st,        St_tpcFieldCage>         ("Geometry/tpc/tpcFieldCage"),
  new ConfigStruct<tpcFieldCageShort_st,   St_tpcFieldCageShort>    ("Geometry/tpc/tpcFieldCageShort"),
  new ConfigStruct<tpcGlobalPosition_st,   St_tpcGlobalPosition>    ("Geometry/tpc/tpcGlobalPosition"),
  new ConfigStruct<Survey_st,              St_Survey>               ("Geometry/tpc/TpcHalfPosition"),
  new ConfigStruct<tpcHVPlanes_st,         St_tpcHVPlanes>          ("Geometry/tpc/tpcHVPlanes"),
  new ConfigStruct<Survey_st,              St_Survey>               ("Geometry/tpc/TpcInnerSectorPosition", "Geometry/tpc/TpcInnerSectorPositionB"),
  new ConfigStruct<Survey_st,              St_Survey>               ("Geometry/tpc/TpcOuterSectorPosition", "Geometry/tpc/TpcOuterSectorPositionB"),
  new ConfigStruct<tpcPadConfig_st,        St_tpcPadConfig>         ("Geometry/tpc/tpcPadConfig"),
  new ConfigStruct<tpcPadPlanes_st,        St_tpcPadPlanes>         ("Geometry/tpc/tpcPadPlanes"),
  new ConfigStruct<Survey_st,              St_Survey>               ("Geometry/tpc/TpcPosition"),
  new ConfigStruct<Survey_st,              St_Survey>               ("Geometry/tpc/TpcSuperSectorPosition", "Geometry/tpc/TpcSuperSectorPositionB"),
  new ConfigStruct<tpcWirePlanes_st,       St_tpcWirePlanes>        ("Geometry/tpc/tpcWirePlanes"),
  new ConfigStruct<MagFactor_st,           St_MagFactor>            ("RunLog/MagFactor"),
  new ConfigStruct<beamInfo_st,            St_beamInfo>             ("RunLog/onl/beamInfo"),
  new ConfigStruct<starClockOnl_st,        St_starClockOnl>         ("RunLog/onl/starClockOnl"),
  new ConfigStruct<starMagOnl_st,          St_starMagOnl>           ("RunLog/onl/starMagOnl"),
  new ConfigStruct<tpcRDOMasks_st,         St_tpcRDOMasks>          ("RunLog/onl/tpcRDOMasks"),
  new ConfigStruct<tss_tsspar_st,          St_tss_tsspar>           ("tpc/tsspars/tsspar")
};

}

#endif
