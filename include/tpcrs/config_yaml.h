#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"
#include "tpcrs/config_type.h"

namespace YAML {
using std::copy;
using std::begin;
using std::end;
using std::array;
using std::string;
using std::vector;
}


namespace YAML {
template<>
struct convert<tpcCalibResolutions> {
  static Node encode(const tpcCalibResolutions& st) {
    Node node;
    
    node["SpaceCharge"] = st.SpaceCharge;
    node["GridLeak"] = st.GridLeak;
    node["comment"] = string(st.comment);
    return node;
  };

  static bool decode(const Node& node, tpcCalibResolutions& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.SpaceCharge = node["SpaceCharge"].as<float>();
    st.GridLeak = node["GridLeak"].as<float>();
    auto comment = node["comment"].as<string>();
    copy(begin(comment), end(comment), st.comment);
    st.comment[comment.size()] ='\0';
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<MDFCorrection> {
  static Node encode(const MDFCorrection& st) {
    Node node;
    
    node["idx"] = static_cast<int>(st.idx);
    node["nrows"] = static_cast<int>(st.nrows);
    node["PolyType"] = static_cast<int>(st.PolyType);
    node["NVariables"] = static_cast<int>(st.NVariables);
    node["NCoefficients"] = static_cast<int>(st.NCoefficients);
    node["Power"] = vector<int>(st.Power, st.Power + 100);
    node["Power"].SetStyle(YAML::EmitterStyle::Flow);
    node["DMean"] = st.DMean;
    node["XMin"] = reinterpret_cast<const array<double, 2>&>( st.XMin );
    node["XMin"].SetStyle(YAML::EmitterStyle::Flow);
    node["XMax"] = reinterpret_cast<const array<double, 2>&>( st.XMax );
    node["XMax"].SetStyle(YAML::EmitterStyle::Flow);
    node["Coefficients"] = reinterpret_cast<const array<double, 50>&>( st.Coefficients );
    node["Coefficients"].SetStyle(YAML::EmitterStyle::Flow);
    node["CoefficientsRMS"] = reinterpret_cast<const array<double, 50>&>( st.CoefficientsRMS );
    node["CoefficientsRMS"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, MDFCorrection& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.idx = static_cast<unsigned char>(node["idx"].as<int>());
    st.nrows = static_cast<unsigned char>(node["nrows"].as<int>());
    st.PolyType = static_cast<unsigned char>(node["PolyType"].as<int>());
    st.NVariables = static_cast<unsigned char>(node["NVariables"].as<int>());
    st.NCoefficients = static_cast<unsigned char>(node["NCoefficients"].as<int>());
    auto Power = node["Power"].as< vector<int> >();
    copy(begin(Power), end(Power), st.Power);
    st.DMean = node["DMean"].as<double>();
    auto XMin = node["XMin"].as<array<double, 2>>();
    copy(begin(XMin), end(XMin), reinterpret_cast<double*>(st.XMin));
    auto XMax = node["XMax"].as<array<double, 2>>();
    copy(begin(XMax), end(XMax), reinterpret_cast<double*>(st.XMax));
    auto Coefficients = node["Coefficients"].as<array<double, 50>>();
    copy(begin(Coefficients), end(Coefficients), reinterpret_cast<double*>(st.Coefficients));
    auto CoefficientsRMS = node["CoefficientsRMS"].as<array<double, 50>>();
    copy(begin(CoefficientsRMS), end(CoefficientsRMS), reinterpret_cast<double*>(st.CoefficientsRMS));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcGlobalPosition> {
  static Node encode(const tpcGlobalPosition& st) {
    Node node;
    
    node["LocalxShift"] = st.LocalxShift;
    node["LocalyShift"] = st.LocalyShift;
    node["LocalzShift"] = st.LocalzShift;
    node["PhiXY"] = st.PhiXY;
    node["PhiXZ"] = st.PhiXZ;
    node["PhiYZ"] = st.PhiYZ;
    node["XX"] = st.XX;
    node["YY"] = st.YY;
    node["ZZ"] = st.ZZ;
    node["PhiXY_geom"] = st.PhiXY_geom;
    node["PhiXZ_geom"] = st.PhiXZ_geom;
    node["PhiYZ_geom"] = st.PhiYZ_geom;
    node["XX_geom"] = st.XX_geom;
    node["YY_geom"] = st.YY_geom;
    node["ZZ_geom"] = st.ZZ_geom;
    return node;
  };

  static bool decode(const Node& node, tpcGlobalPosition& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.LocalxShift = node["LocalxShift"].as<float>();
    st.LocalyShift = node["LocalyShift"].as<float>();
    st.LocalzShift = node["LocalzShift"].as<float>();
    st.PhiXY = node["PhiXY"].as<float>();
    st.PhiXZ = node["PhiXZ"].as<float>();
    st.PhiYZ = node["PhiYZ"].as<float>();
    st.XX = node["XX"].as<float>();
    st.YY = node["YY"].as<float>();
    st.ZZ = node["ZZ"].as<float>();
    st.PhiXY_geom = node["PhiXY_geom"].as<float>();
    st.PhiXZ_geom = node["PhiXZ_geom"].as<float>();
    st.PhiYZ_geom = node["PhiYZ_geom"].as<float>();
    st.XX_geom = node["XX_geom"].as<float>();
    st.YY_geom = node["YY_geom"].as<float>();
    st.ZZ_geom = node["ZZ_geom"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcRDOMasks> {
  static Node encode(const tpcRDOMasks& st) {
    Node node;
    
    node["runNumber"] = st.runNumber;
    node["sector"] = st.sector;
    node["mask"] = st.mask;
    return node;
  };

  static bool decode(const Node& node, tpcRDOMasks& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.runNumber = node["runNumber"].as<unsigned int>();
    st.sector = node["sector"].as<unsigned int>();
    st.mask = node["mask"].as<unsigned int>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tss_tsspar> {
  static Node encode(const tss_tsspar& st) {
    Node node;
    
    node["fileout"] = string(st.fileout);
    node["dynam"] = st.dynam;
    node["format"] = st.format;
    node["max_itime"] = st.max_itime;
    node["max_pads"] = st.max_pads;
    node["max_row"] = st.max_row;
    node["max_sect"] = st.max_sect;
    node["min_itime"] = st.min_itime;
    node["min_pads"] = st.min_pads;
    node["min_row"] = st.min_row;
    node["min_sect"] = st.min_sect;
    node["mode"] = st.mode;
    node["nele_laser"] = st.nele_laser;
    node["ngain"] = st.ngain;
    node["nseg"] = st.nseg;
    node["ntime"] = st.ntime;
    node["printout"] = st.printout;
    node["tpc_half"] = st.tpc_half;
    node["reset"] = st.reset;
    node["ave_ion_pot"] = st.ave_ion_pot;
    node["bfield"] = st.bfield;
    node["c_test"] = st.c_test;
    node["diff_long"] = st.diff_long;
    node["diff_trans"] = st.diff_trans;
    node["gain_in"] = st.gain_in;
    node["gain_out"] = st.gain_out;
    node["prf_in"] = st.prf_in;
    node["prf_out"] = st.prf_out;
    node["sca_rms"] = st.sca_rms;
    node["scale"] = st.scale;
    node["step_size"] = st.step_size;
    node["tau"] = st.tau;
    node["threshold"] = st.threshold;
    node["time_offset"] = st.time_offset;
    node["v_test"] = st.v_test;
    node["white_rms"] = st.white_rms;
    node["wire_coupling_in"] = st.wire_coupling_in;
    node["wire_coupling_out"] = st.wire_coupling_out;
    node["x_laser"] = st.x_laser;
    node["y_laser"] = st.y_laser;
    node["z_laser"] = st.z_laser;
    return node;
  };

  static bool decode(const Node& node, tss_tsspar& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    auto fileout = node["fileout"].as<string>();
    copy(begin(fileout), end(fileout), st.fileout);
    st.fileout[fileout.size()] ='\0';
    st.dynam = node["dynam"].as<int>();
    st.format = node["format"].as<int>();
    st.max_itime = node["max_itime"].as<int>();
    st.max_pads = node["max_pads"].as<int>();
    st.max_row = node["max_row"].as<int>();
    st.max_sect = node["max_sect"].as<int>();
    st.min_itime = node["min_itime"].as<int>();
    st.min_pads = node["min_pads"].as<int>();
    st.min_row = node["min_row"].as<int>();
    st.min_sect = node["min_sect"].as<int>();
    st.mode = node["mode"].as<int>();
    st.nele_laser = node["nele_laser"].as<int>();
    st.ngain = node["ngain"].as<int>();
    st.nseg = node["nseg"].as<int>();
    st.ntime = node["ntime"].as<int>();
    st.printout = node["printout"].as<int>();
    st.tpc_half = node["tpc_half"].as<int>();
    st.reset = node["reset"].as<int>();
    st.ave_ion_pot = node["ave_ion_pot"].as<float>();
    st.bfield = node["bfield"].as<float>();
    st.c_test = node["c_test"].as<float>();
    st.diff_long = node["diff_long"].as<float>();
    st.diff_trans = node["diff_trans"].as<float>();
    st.gain_in = node["gain_in"].as<float>();
    st.gain_out = node["gain_out"].as<float>();
    st.prf_in = node["prf_in"].as<float>();
    st.prf_out = node["prf_out"].as<float>();
    st.sca_rms = node["sca_rms"].as<float>();
    st.scale = node["scale"].as<float>();
    st.step_size = node["step_size"].as<float>();
    st.tau = node["tau"].as<float>();
    st.threshold = node["threshold"].as<float>();
    st.time_offset = node["time_offset"].as<float>();
    st.v_test = node["v_test"].as<float>();
    st.white_rms = node["white_rms"].as<float>();
    st.wire_coupling_in = node["wire_coupling_in"].as<float>();
    st.wire_coupling_out = node["wire_coupling_out"].as<float>();
    st.x_laser = node["x_laser"].as<float>();
    st.y_laser = node["y_laser"].as<float>();
    st.z_laser = node["z_laser"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcOmegaTau> {
  static Node encode(const tpcOmegaTau& st) {
    Node node;
    
    node["tensorV1"] = st.tensorV1;
    node["tensorV2"] = st.tensorV2;
    node["distortionCorrectionsMode"] = st.distortionCorrectionsMode;
    return node;
  };

  static bool decode(const Node& node, tpcOmegaTau& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.tensorV1 = node["tensorV1"].as<float>();
    st.tensorV2 = node["tensorV2"].as<float>();
    st.distortionCorrectionsMode = node["distortionCorrectionsMode"].as<unsigned short>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcFieldCageShort> {
  static Node encode(const tpcFieldCageShort& st) {
    Node node;
    
    node["side"] = st.side;
    node["cage"] = st.cage;
    node["ring"] = st.ring;
    node["resistor"] = st.resistor;
    node["MissingResistance"] = st.MissingResistance;
    return node;
  };

  static bool decode(const Node& node, tpcFieldCageShort& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.side = node["side"].as<float>();
    st.cage = node["cage"].as<float>();
    st.ring = node["ring"].as<float>();
    st.resistor = node["resistor"].as<float>();
    st.MissingResistance = node["MissingResistance"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<g2t_tpc_hit> {
  static Node encode(const g2t_tpc_hit& st) {
    Node node;
    
    node["id"] = st.id;
    node["next_tr_hit_p"] = st.next_tr_hit_p;
    node["track_p"] = st.track_p;
    node["volume_id"] = st.volume_id;
    node["de"] = st.de;
    node["ds"] = st.ds;
    node["p"] = reinterpret_cast<const array<float, 3>&>( st.p );
    node["p"].SetStyle(YAML::EmitterStyle::Flow);
    node["tof"] = st.tof;
    node["x"] = reinterpret_cast<const array<float, 3>&>( st.x );
    node["x"].SetStyle(YAML::EmitterStyle::Flow);
    node["lgam"] = st.lgam;
    node["length"] = st.length;
    node["adc"] = st.adc;
    node["pad"] = st.pad;
    node["timebucket"] = st.timebucket;
    node["np"] = st.np;
    return node;
  };

  static bool decode(const Node& node, g2t_tpc_hit& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.id = node["id"].as<int>();
    st.next_tr_hit_p = node["next_tr_hit_p"].as<int>();
    st.track_p = node["track_p"].as<int>();
    st.volume_id = node["volume_id"].as<int>();
    st.de = node["de"].as<float>();
    st.ds = node["ds"].as<float>();
    auto p = node["p"].as<array<float, 3>>();
    copy(begin(p), end(p), reinterpret_cast<float*>(st.p));
    st.tof = node["tof"].as<float>();
    auto x = node["x"].as<array<float, 3>>();
    copy(begin(x), end(x), reinterpret_cast<float*>(st.x));
    st.lgam = node["lgam"].as<float>();
    st.length = node["length"].as<float>();
    st.adc = node["adc"].as<float>();
    st.pad = node["pad"].as<float>();
    st.timebucket = node["timebucket"].as<float>();
    st.np = node["np"].as<int>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<TpcSecRowCor> {
  static Node encode(const TpcSecRowCor& st) {
    Node node;
    
    node["GainScale"] = reinterpret_cast<const array<float, 100>&>( st.GainScale );
    node["GainScale"].SetStyle(YAML::EmitterStyle::Flow);
    node["GainRms"] = reinterpret_cast<const array<float, 100>&>( st.GainRms );
    node["GainRms"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, TpcSecRowCor& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    auto GainScale = node["GainScale"].as<array<float, 100>>();
    copy(begin(GainScale), end(GainScale), reinterpret_cast<float*>(st.GainScale));
    auto GainRms = node["GainRms"].as<array<float, 100>>();
    copy(begin(GainRms), end(GainRms), reinterpret_cast<float*>(st.GainRms));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcRDOMap> {
  static Node encode(const tpcRDOMap& st) {
    Node node;
    
    node["idx"] = st.idx;
    node["nrows"] = st.nrows;
    node["row"] = static_cast<int>(st.row);
    node["padMin"] = static_cast<int>(st.padMin);
    node["padMax"] = static_cast<int>(st.padMax);
    node["rdo"] = static_cast<int>(st.rdo);
    return node;
  };

  static bool decode(const Node& node, tpcRDOMap& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.idx = node["idx"].as<int>();
    st.nrows = node["nrows"].as<int>();
    st.row = static_cast<unsigned char>(node["row"].as<int>());
    st.padMin = static_cast<unsigned char>(node["padMin"].as<int>());
    st.padMax = static_cast<unsigned char>(node["padMax"].as<int>());
    st.rdo = static_cast<unsigned char>(node["rdo"].as<int>());
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcAnodeHVavg> {
  static Node encode(const tpcAnodeHVavg& st) {
    Node node;
    
    node["sector"] = st.sector;
    node["socket"] = st.socket;
    node["voltage"] = st.voltage;
    node["rms"] = st.rms;
    node["numentries"] = st.numentries;
    node["numoutliers"] = st.numoutliers;
    return node;
  };

  static bool decode(const Node& node, tpcAnodeHVavg& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.sector = node["sector"].as<unsigned short>();
    st.socket = node["socket"].as<unsigned short>();
    st.voltage = node["voltage"].as<float>();
    st.rms = node["rms"].as<float>();
    st.numentries = node["numentries"].as<int>();
    st.numoutliers = node["numoutliers"].as<int>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcFieldCage> {
  static Node encode(const tpcFieldCage& st) {
    Node node;
    
    node["innerFieldCageShift"] = st.innerFieldCageShift;
    node["eastClockError"] = st.eastClockError;
    node["westClockError"] = st.westClockError;
    return node;
  };

  static bool decode(const Node& node, tpcFieldCage& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.innerFieldCageShift = node["innerFieldCageShift"].as<float>();
    st.eastClockError = node["eastClockError"].as<float>();
    st.westClockError = node["westClockError"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<trigDetSums> {
  static Node encode(const trigDetSums& st) {
    Node node;
    
    node["runNumber"] = st.runNumber;
    node["timeOffset"] = st.timeOffset;
    node["ctbWest"] = st.ctbWest;
    node["ctbEast"] = st.ctbEast;
    node["ctbTOFp"] = st.ctbTOFp;
    node["tofp"] = st.tofp;
    node["zdcWest"] = st.zdcWest;
    node["zdcEast"] = st.zdcEast;
    node["zdcX"] = st.zdcX;
    node["mult"] = st.mult;
    node["L0"] = st.L0;
    node["bbcX"] = st.bbcX;
    node["bbcXctbTOFp"] = st.bbcXctbTOFp;
    node["bbcWest"] = st.bbcWest;
    node["bbcEast"] = st.bbcEast;
    node["bbcYellowBkg"] = st.bbcYellowBkg;
    node["bbcBlueBkg"] = st.bbcBlueBkg;
    node["pvpdWest"] = st.pvpdWest;
    node["pvpdEast"] = st.pvpdEast;
    return node;
  };

  static bool decode(const Node& node, trigDetSums& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.runNumber = node["runNumber"].as<unsigned int>();
    st.timeOffset = node["timeOffset"].as<unsigned int>();
    st.ctbWest = node["ctbWest"].as<double>();
    st.ctbEast = node["ctbEast"].as<double>();
    st.ctbTOFp = node["ctbTOFp"].as<double>();
    st.tofp = node["tofp"].as<double>();
    st.zdcWest = node["zdcWest"].as<double>();
    st.zdcEast = node["zdcEast"].as<double>();
    st.zdcX = node["zdcX"].as<double>();
    st.mult = node["mult"].as<double>();
    st.L0 = node["L0"].as<double>();
    st.bbcX = node["bbcX"].as<double>();
    st.bbcXctbTOFp = node["bbcXctbTOFp"].as<double>();
    st.bbcWest = node["bbcWest"].as<double>();
    st.bbcEast = node["bbcEast"].as<double>();
    st.bbcYellowBkg = node["bbcYellowBkg"].as<double>();
    st.bbcBlueBkg = node["bbcBlueBkg"].as<double>();
    st.pvpdWest = node["pvpdWest"].as<double>();
    st.pvpdEast = node["pvpdEast"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcPadPlanes> {
  static Node encode(const tpcPadPlanes& st) {
    Node node;
    
    node["padRows"] = st.padRows;
    node["innerPadRows"] = st.innerPadRows;
    node["innerPadRows48"] = st.innerPadRows48;
    node["innerPadRows52"] = st.innerPadRows52;
    node["outerPadRows"] = st.outerPadRows;
    node["superInnerPadRows"] = st.superInnerPadRows;
    node["superOuterPadRows"] = st.superOuterPadRows;
    node["innerSectorPadWidth"] = st.innerSectorPadWidth;
    node["innerSectorPadLength"] = st.innerSectorPadLength;
    node["innerSectorPadPitch"] = st.innerSectorPadPitch;
    node["innerSectorRowPitch1"] = st.innerSectorRowPitch1;
    node["innerSectorRowPitch2"] = st.innerSectorRowPitch2;
    node["firstPadRow"] = st.firstPadRow;
    node["firstOuterSectorPadRow"] = st.firstOuterSectorPadRow;
    node["lastOuterSectorPadRow"] = st.lastOuterSectorPadRow;
    node["firstRowWidth"] = st.firstRowWidth;
    node["lastRowWidth"] = st.lastRowWidth;
    node["outerSectorPadWidth"] = st.outerSectorPadWidth;
    node["outerSectorPadLength"] = st.outerSectorPadLength;
    node["outerSectorPadPitch"] = st.outerSectorPadPitch;
    node["outerSectorRowPitch"] = st.outerSectorRowPitch;
    node["outerSectorLength"] = st.outerSectorLength;
    node["ioSectorSeparation"] = st.ioSectorSeparation;
    node["innerSectorEdge"] = st.innerSectorEdge;
    node["outerSectorEdge"] = st.outerSectorEdge;
    node["innerSectorPadPlaneZ"] = st.innerSectorPadPlaneZ;
    node["outerSectorPadPlaneZ"] = st.outerSectorPadPlaneZ;
    node["innerPadsPerRow"] = reinterpret_cast<const array<int, 13>&>( st.innerPadsPerRow );
    node["innerPadsPerRow"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerPadsPerRow"] = reinterpret_cast<const array<int, 32>&>( st.outerPadsPerRow );
    node["outerPadsPerRow"].SetStyle(YAML::EmitterStyle::Flow);
    node["innerRowRadii"] = reinterpret_cast<const array<double, 13>&>( st.innerRowRadii );
    node["innerRowRadii"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerRowRadii"] = reinterpret_cast<const array<double, 32>&>( st.outerRowRadii );
    node["outerRowRadii"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcPadPlanes& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.padRows = node["padRows"].as<int>();
    st.innerPadRows = node["innerPadRows"].as<int>();
    st.innerPadRows48 = node["innerPadRows48"].as<int>();
    st.innerPadRows52 = node["innerPadRows52"].as<int>();
    st.outerPadRows = node["outerPadRows"].as<int>();
    st.superInnerPadRows = node["superInnerPadRows"].as<int>();
    st.superOuterPadRows = node["superOuterPadRows"].as<int>();
    st.innerSectorPadWidth = node["innerSectorPadWidth"].as<double>();
    st.innerSectorPadLength = node["innerSectorPadLength"].as<double>();
    st.innerSectorPadPitch = node["innerSectorPadPitch"].as<double>();
    st.innerSectorRowPitch1 = node["innerSectorRowPitch1"].as<double>();
    st.innerSectorRowPitch2 = node["innerSectorRowPitch2"].as<double>();
    st.firstPadRow = node["firstPadRow"].as<double>();
    st.firstOuterSectorPadRow = node["firstOuterSectorPadRow"].as<double>();
    st.lastOuterSectorPadRow = node["lastOuterSectorPadRow"].as<double>();
    st.firstRowWidth = node["firstRowWidth"].as<double>();
    st.lastRowWidth = node["lastRowWidth"].as<double>();
    st.outerSectorPadWidth = node["outerSectorPadWidth"].as<double>();
    st.outerSectorPadLength = node["outerSectorPadLength"].as<double>();
    st.outerSectorPadPitch = node["outerSectorPadPitch"].as<double>();
    st.outerSectorRowPitch = node["outerSectorRowPitch"].as<double>();
    st.outerSectorLength = node["outerSectorLength"].as<double>();
    st.ioSectorSeparation = node["ioSectorSeparation"].as<double>();
    st.innerSectorEdge = node["innerSectorEdge"].as<double>();
    st.outerSectorEdge = node["outerSectorEdge"].as<double>();
    st.innerSectorPadPlaneZ = node["innerSectorPadPlaneZ"].as<double>();
    st.outerSectorPadPlaneZ = node["outerSectorPadPlaneZ"].as<double>();
    auto innerPadsPerRow = node["innerPadsPerRow"].as<array<int, 13>>();
    copy(begin(innerPadsPerRow), end(innerPadsPerRow), reinterpret_cast<int*>(st.innerPadsPerRow));
    auto outerPadsPerRow = node["outerPadsPerRow"].as<array<int, 32>>();
    copy(begin(outerPadsPerRow), end(outerPadsPerRow), reinterpret_cast<int*>(st.outerPadsPerRow));
    auto innerRowRadii = node["innerRowRadii"].as<array<double, 13>>();
    copy(begin(innerRowRadii), end(innerRowRadii), reinterpret_cast<double*>(st.innerRowRadii));
    auto outerRowRadii = node["outerRowRadii"].as<array<double, 32>>();
    copy(begin(outerRowRadii), end(outerRowRadii), reinterpret_cast<double*>(st.outerRowRadii));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcDriftVelocity> {
  static Node encode(const tpcDriftVelocity& st) {
    Node node;
    
    node["laserDriftVelocityEast"] = st.laserDriftVelocityEast;
    node["laserDriftVelocityWest"] = st.laserDriftVelocityWest;
    node["cathodeDriftVelocityEast"] = st.cathodeDriftVelocityEast;
    node["cathodeDriftVelocityWest"] = st.cathodeDriftVelocityWest;
    return node;
  };

  static bool decode(const Node& node, tpcDriftVelocity& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.laserDriftVelocityEast = node["laserDriftVelocityEast"].as<float>();
    st.laserDriftVelocityWest = node["laserDriftVelocityWest"].as<float>();
    st.cathodeDriftVelocityEast = node["cathodeDriftVelocityEast"].as<float>();
    st.cathodeDriftVelocityWest = node["cathodeDriftVelocityWest"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<starClockOnl> {
  static Node encode(const starClockOnl& st) {
    Node node;
    
    node["runNumber"] = st.runNumber;
    node["time"] = st.time;
    node["frequency"] = st.frequency;
    return node;
  };

  static bool decode(const Node& node, starClockOnl& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.runNumber = node["runNumber"].as<unsigned int>();
    st.time = node["time"].as<unsigned int>();
    st.frequency = node["frequency"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcWirePlanes> {
  static Node encode(const tpcWirePlanes& st) {
    Node node;
    
    node["anodeWireRadius"] = st.anodeWireRadius;
    node["frischGridWireRadius"] = st.frischGridWireRadius;
    node["gatingGridWireRadius"] = st.gatingGridWireRadius;
    node["anodeWirePitch"] = st.anodeWirePitch;
    node["frischGridWirePitch"] = st.frischGridWirePitch;
    node["gatingGridWirePitch"] = st.gatingGridWirePitch;
    node["innerSectorAnodeWirePadSep"] = st.innerSectorAnodeWirePadSep;
    node["innerSectorFrischGridPadSep"] = st.innerSectorFrischGridPadSep;
    node["innerSectorGatingGridPadSep"] = st.innerSectorGatingGridPadSep;
    node["outerSectorAnodeWirePadSep"] = st.outerSectorAnodeWirePadSep;
    node["outerSectorFrischGridPadSep"] = st.outerSectorFrischGridPadSep;
    node["outerSectorGatingGridPadSep"] = st.outerSectorGatingGridPadSep;
    node["numInnerSectorAnodeWires"] = st.numInnerSectorAnodeWires;
    node["numInnerSectorFrischGridWires"] = st.numInnerSectorFrischGridWires;
    node["numInnerSectorGatingGridWires"] = st.numInnerSectorGatingGridWires;
    node["firstInnerSectorAnodeWire"] = st.firstInnerSectorAnodeWire;
    node["firstInnerSectorFrischGridWire"] = st.firstInnerSectorFrischGridWire;
    node["firstInnerSectorGatingGridWire"] = st.firstInnerSectorGatingGridWire;
    node["lastInnerSectorAnodeWire"] = st.lastInnerSectorAnodeWire;
    node["numOuterSectorAnodeWires"] = st.numOuterSectorAnodeWires;
    node["numOuterSectorFrischGridWires"] = st.numOuterSectorFrischGridWires;
    node["numOuterSectorGatingGridWires"] = st.numOuterSectorGatingGridWires;
    node["firstOuterSectorAnodeWire"] = st.firstOuterSectorAnodeWire;
    node["firstOuterSectorFrischGridWire"] = st.firstOuterSectorFrischGridWire;
    node["firstOuterSectorGatingGridWire"] = st.firstOuterSectorGatingGridWire;
    node["lastOuterSectorAnodeWire"] = st.lastOuterSectorAnodeWire;
    return node;
  };

  static bool decode(const Node& node, tpcWirePlanes& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.anodeWireRadius = node["anodeWireRadius"].as<double>();
    st.frischGridWireRadius = node["frischGridWireRadius"].as<double>();
    st.gatingGridWireRadius = node["gatingGridWireRadius"].as<double>();
    st.anodeWirePitch = node["anodeWirePitch"].as<double>();
    st.frischGridWirePitch = node["frischGridWirePitch"].as<double>();
    st.gatingGridWirePitch = node["gatingGridWirePitch"].as<double>();
    st.innerSectorAnodeWirePadSep = node["innerSectorAnodeWirePadSep"].as<double>();
    st.innerSectorFrischGridPadSep = node["innerSectorFrischGridPadSep"].as<double>();
    st.innerSectorGatingGridPadSep = node["innerSectorGatingGridPadSep"].as<double>();
    st.outerSectorAnodeWirePadSep = node["outerSectorAnodeWirePadSep"].as<double>();
    st.outerSectorFrischGridPadSep = node["outerSectorFrischGridPadSep"].as<double>();
    st.outerSectorGatingGridPadSep = node["outerSectorGatingGridPadSep"].as<double>();
    st.numInnerSectorAnodeWires = node["numInnerSectorAnodeWires"].as<int>();
    st.numInnerSectorFrischGridWires = node["numInnerSectorFrischGridWires"].as<int>();
    st.numInnerSectorGatingGridWires = node["numInnerSectorGatingGridWires"].as<int>();
    st.firstInnerSectorAnodeWire = node["firstInnerSectorAnodeWire"].as<double>();
    st.firstInnerSectorFrischGridWire = node["firstInnerSectorFrischGridWire"].as<double>();
    st.firstInnerSectorGatingGridWire = node["firstInnerSectorGatingGridWire"].as<double>();
    st.lastInnerSectorAnodeWire = node["lastInnerSectorAnodeWire"].as<double>();
    st.numOuterSectorAnodeWires = node["numOuterSectorAnodeWires"].as<int>();
    st.numOuterSectorFrischGridWires = node["numOuterSectorFrischGridWires"].as<int>();
    st.numOuterSectorGatingGridWires = node["numOuterSectorGatingGridWires"].as<int>();
    st.firstOuterSectorAnodeWire = node["firstOuterSectorAnodeWire"].as<double>();
    st.firstOuterSectorFrischGridWire = node["firstOuterSectorFrischGridWire"].as<double>();
    st.firstOuterSectorGatingGridWire = node["firstOuterSectorGatingGridWire"].as<double>();
    st.lastOuterSectorAnodeWire = node["lastOuterSectorAnodeWire"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcSectorT0offset> {
  static Node encode(const tpcSectorT0offset& st) {
    Node node;
    
    node["t0"] = reinterpret_cast<const array<float, 48>&>( st.t0 );
    node["t0"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcSectorT0offset& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    auto t0 = node["t0"].as<array<float, 48>>();
    copy(begin(t0), end(t0), reinterpret_cast<float*>(st.t0));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcHighVoltages> {
  static Node encode(const tpcHighVoltages& st) {
    Node node;
    
    node["cathode"] = st.cathode;
    node["gatedGridRef"] = st.gatedGridRef;
    node["gridLeakWallTip"] = reinterpret_cast<const array<float, 24>&>( st.gridLeakWallTip );
    node["gridLeakWallTip"].SetStyle(YAML::EmitterStyle::Flow);
    node["gridLeakWallSide"] = reinterpret_cast<const array<float, 24>&>( st.gridLeakWallSide );
    node["gridLeakWallSide"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcHighVoltages& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.cathode = node["cathode"].as<float>();
    st.gatedGridRef = node["gatedGridRef"].as<float>();
    auto gridLeakWallTip = node["gridLeakWallTip"].as<array<float, 24>>();
    copy(begin(gridLeakWallTip), end(gridLeakWallTip), reinterpret_cast<float*>(st.gridLeakWallTip));
    auto gridLeakWallSide = node["gridLeakWallSide"].as<array<float, 24>>();
    copy(begin(gridLeakWallSide), end(gridLeakWallSide), reinterpret_cast<float*>(st.gridLeakWallSide));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcCorrection> {
  static Node encode(const tpcCorrection& st) {
    Node node;
    
    node["type"] = st.type;
    node["idx"] = st.idx;
    node["nrows"] = st.nrows;
    node["npar"] = st.npar;
    node["OffSet"] = st.OffSet;
    node["min"] = st.min;
    node["max"] = st.max;
    node["a"] = reinterpret_cast<const array<double, 10>&>( st.a );
    node["a"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcCorrection& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.type = node["type"].as<int>();
    st.idx = node["idx"].as<int>();
    st.nrows = node["nrows"].as<int>();
    st.npar = node["npar"].as<int>();
    st.OffSet = node["OffSet"].as<double>();
    st.min = node["min"].as<double>();
    st.max = node["max"].as<double>();
    auto a = node["a"].as<array<double, 10>>();
    copy(begin(a), end(a), reinterpret_cast<double*>(st.a));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<g2t_track> {
  static Node encode(const g2t_track& st) {
    Node node;
    
    node["eg_label"] = st.eg_label;
    node["eg_pid"] = st.eg_pid;
    node["ge_pid"] = st.ge_pid;
    node["hit_ctb_p"] = st.hit_ctb_p;
    node["hit_eem_p"] = st.hit_eem_p;
    node["hit_emc_p"] = st.hit_emc_p;
    node["hit_esm_p"] = st.hit_esm_p;
    node["hit_ftp_p"] = st.hit_ftp_p;
    node["hit_gem_p"] = st.hit_gem_p;
    node["hit_hpd_p"] = st.hit_hpd_p;
    node["hit_ist_p"] = st.hit_ist_p;
    node["hit_igt_p"] = st.hit_igt_p;
    node["hit_fst_p"] = st.hit_fst_p;
    node["hit_fgt_p"] = st.hit_fgt_p;
    node["hit_fpd_p"] = st.hit_fpd_p;
    node["hit_fsc_p"] = st.hit_fsc_p;
    node["hit_mtd_p"] = st.hit_mtd_p;
    node["hit_mwc_p"] = st.hit_mwc_p;
    node["hit_pgc_p"] = st.hit_pgc_p;
    node["hit_pmd_p"] = st.hit_pmd_p;
    node["hit_smd_p"] = st.hit_smd_p;
    node["hit_ssd_p"] = st.hit_ssd_p;
    node["hit_svt_p"] = st.hit_svt_p;
    node["hit_pix_p"] = st.hit_pix_p;
    node["hit_tof_p"] = st.hit_tof_p;
    node["hit_tpc_p"] = st.hit_tpc_p;
    node["hit_vpd_p"] = st.hit_vpd_p;
    node["hit_etr_p"] = st.hit_etr_p;
    node["hit_hca_p"] = st.hit_hca_p;
    node["hit_fts_p"] = st.hit_fts_p;
    node["hit_eto_p"] = st.hit_eto_p;
    node["id"] = st.id;
    node["is_shower"] = st.is_shower;
    node["itrmd_vertex_p"] = st.itrmd_vertex_p;
    node["n_ctb_hit"] = st.n_ctb_hit;
    node["n_eem_hit"] = st.n_eem_hit;
    node["n_emc_hit"] = st.n_emc_hit;
    node["n_esm_hit"] = st.n_esm_hit;
    node["n_ftp_hit"] = st.n_ftp_hit;
    node["n_gem_hit"] = st.n_gem_hit;
    node["n_hpd_hit"] = st.n_hpd_hit;
    node["n_ist_hit"] = st.n_ist_hit;
    node["n_igt_hit"] = st.n_igt_hit;
    node["n_fst_hit"] = st.n_fst_hit;
    node["n_fgt_hit"] = st.n_fgt_hit;
    node["n_fpd_hit"] = st.n_fpd_hit;
    node["n_fsc_hit"] = st.n_fsc_hit;
    node["n_mtd_hit"] = st.n_mtd_hit;
    node["n_mwc_hit"] = st.n_mwc_hit;
    node["n_pgc_hit"] = st.n_pgc_hit;
    node["n_pmd_hit"] = st.n_pmd_hit;
    node["n_smd_hit"] = st.n_smd_hit;
    node["n_ssd_hit"] = st.n_ssd_hit;
    node["n_svt_hit"] = st.n_svt_hit;
    node["n_pix_hit"] = st.n_pix_hit;
    node["n_tof_hit"] = st.n_tof_hit;
    node["n_tpc_hit"] = st.n_tpc_hit;
    node["n_vpd_hit"] = st.n_vpd_hit;
    node["n_etr_hit"] = st.n_etr_hit;
    node["n_hca_hit"] = st.n_hca_hit;
    node["n_fts_hit"] = st.n_fts_hit;
    node["n_eto_hit"] = st.n_eto_hit;
    node["n_stg_hit"] = st.n_stg_hit;
    node["n_wca_hit"] = st.n_wca_hit;
    node["next_parent_p"] = st.next_parent_p;
    node["next_vtx_trk_p"] = st.next_vtx_trk_p;
    node["start_vertex_p"] = st.start_vertex_p;
    node["stop_vertex_p"] = st.stop_vertex_p;
    node["charge"] = st.charge;
    node["e"] = st.e;
    node["eta"] = st.eta;
    node["p"] = reinterpret_cast<const array<float, 3>&>( st.p );
    node["p"].SetStyle(YAML::EmitterStyle::Flow);
    node["pt"] = st.pt;
    node["ptot"] = st.ptot;
    node["rapidity"] = st.rapidity;
    return node;
  };

  static bool decode(const Node& node, g2t_track& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.eg_label = node["eg_label"].as<int>();
    st.eg_pid = node["eg_pid"].as<int>();
    st.ge_pid = node["ge_pid"].as<int>();
    st.hit_ctb_p = node["hit_ctb_p"].as<int>();
    st.hit_eem_p = node["hit_eem_p"].as<int>();
    st.hit_emc_p = node["hit_emc_p"].as<int>();
    st.hit_esm_p = node["hit_esm_p"].as<int>();
    st.hit_ftp_p = node["hit_ftp_p"].as<int>();
    st.hit_gem_p = node["hit_gem_p"].as<int>();
    st.hit_hpd_p = node["hit_hpd_p"].as<int>();
    st.hit_ist_p = node["hit_ist_p"].as<int>();
    st.hit_igt_p = node["hit_igt_p"].as<int>();
    st.hit_fst_p = node["hit_fst_p"].as<int>();
    st.hit_fgt_p = node["hit_fgt_p"].as<int>();
    st.hit_fpd_p = node["hit_fpd_p"].as<int>();
    st.hit_fsc_p = node["hit_fsc_p"].as<int>();
    st.hit_mtd_p = node["hit_mtd_p"].as<int>();
    st.hit_mwc_p = node["hit_mwc_p"].as<int>();
    st.hit_pgc_p = node["hit_pgc_p"].as<int>();
    st.hit_pmd_p = node["hit_pmd_p"].as<int>();
    st.hit_smd_p = node["hit_smd_p"].as<int>();
    st.hit_ssd_p = node["hit_ssd_p"].as<int>();
    st.hit_svt_p = node["hit_svt_p"].as<int>();
    st.hit_pix_p = node["hit_pix_p"].as<int>();
    st.hit_tof_p = node["hit_tof_p"].as<int>();
    st.hit_tpc_p = node["hit_tpc_p"].as<int>();
    st.hit_vpd_p = node["hit_vpd_p"].as<int>();
    st.hit_etr_p = node["hit_etr_p"].as<int>();
    st.hit_hca_p = node["hit_hca_p"].as<int>();
    st.hit_fts_p = node["hit_fts_p"].as<int>();
    st.hit_eto_p = node["hit_eto_p"].as<int>();
    st.id = node["id"].as<int>();
    st.is_shower = node["is_shower"].as<int>();
    st.itrmd_vertex_p = node["itrmd_vertex_p"].as<int>();
    st.n_ctb_hit = node["n_ctb_hit"].as<int>();
    st.n_eem_hit = node["n_eem_hit"].as<int>();
    st.n_emc_hit = node["n_emc_hit"].as<int>();
    st.n_esm_hit = node["n_esm_hit"].as<int>();
    st.n_ftp_hit = node["n_ftp_hit"].as<int>();
    st.n_gem_hit = node["n_gem_hit"].as<int>();
    st.n_hpd_hit = node["n_hpd_hit"].as<int>();
    st.n_ist_hit = node["n_ist_hit"].as<int>();
    st.n_igt_hit = node["n_igt_hit"].as<int>();
    st.n_fst_hit = node["n_fst_hit"].as<int>();
    st.n_fgt_hit = node["n_fgt_hit"].as<int>();
    st.n_fpd_hit = node["n_fpd_hit"].as<int>();
    st.n_fsc_hit = node["n_fsc_hit"].as<int>();
    st.n_mtd_hit = node["n_mtd_hit"].as<int>();
    st.n_mwc_hit = node["n_mwc_hit"].as<int>();
    st.n_pgc_hit = node["n_pgc_hit"].as<int>();
    st.n_pmd_hit = node["n_pmd_hit"].as<int>();
    st.n_smd_hit = node["n_smd_hit"].as<int>();
    st.n_ssd_hit = node["n_ssd_hit"].as<int>();
    st.n_svt_hit = node["n_svt_hit"].as<int>();
    st.n_pix_hit = node["n_pix_hit"].as<int>();
    st.n_tof_hit = node["n_tof_hit"].as<int>();
    st.n_tpc_hit = node["n_tpc_hit"].as<int>();
    st.n_vpd_hit = node["n_vpd_hit"].as<int>();
    st.n_etr_hit = node["n_etr_hit"].as<int>();
    st.n_hca_hit = node["n_hca_hit"].as<int>();
    st.n_fts_hit = node["n_fts_hit"].as<int>();
    st.n_eto_hit = node["n_eto_hit"].as<int>();
    st.n_stg_hit = node["n_stg_hit"].as<int>();
    st.n_wca_hit = node["n_wca_hit"].as<int>();
    st.next_parent_p = node["next_parent_p"].as<int>();
    st.next_vtx_trk_p = node["next_vtx_trk_p"].as<int>();
    st.start_vertex_p = node["start_vertex_p"].as<int>();
    st.stop_vertex_p = node["stop_vertex_p"].as<int>();
    st.charge = node["charge"].as<float>();
    st.e = node["e"].as<float>();
    st.eta = node["eta"].as<float>();
    auto p = node["p"].as<array<float, 3>>();
    copy(begin(p), end(p), reinterpret_cast<float*>(st.p));
    st.pt = node["pt"].as<float>();
    st.ptot = node["ptot"].as<float>();
    st.rapidity = node["rapidity"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<TpcEffectivedX> {
  static Node encode(const TpcEffectivedX& st) {
    Node node;
    
    node["scaleInner"] = st.scaleInner;
    node["scaleOuter"] = st.scaleOuter;
    return node;
  };

  static bool decode(const Node& node, TpcEffectivedX& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.scaleInner = node["scaleInner"].as<float>();
    st.scaleOuter = node["scaleOuter"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcSCGL> {
  static Node encode(const tpcSCGL& st) {
    Node node;
    
    node["SC"] = reinterpret_cast<const array<float, 8>&>( st.SC );
    node["SC"].SetStyle(YAML::EmitterStyle::Flow);
    node["SCoffset"] = reinterpret_cast<const array<float, 8>&>( st.SCoffset );
    node["SCoffset"].SetStyle(YAML::EmitterStyle::Flow);
    node["SCexponent"] = reinterpret_cast<const array<float, 8>&>( st.SCexponent );
    node["SCexponent"].SetStyle(YAML::EmitterStyle::Flow);
    node["SCscaler"] = reinterpret_cast<const array<float, 8>&>( st.SCscaler );
    node["SCscaler"].SetStyle(YAML::EmitterStyle::Flow);
    node["GL"] = reinterpret_cast<const array<float, 24>&>( st.GL );
    node["GL"].SetStyle(YAML::EmitterStyle::Flow);
    node["GLoffset"] = reinterpret_cast<const array<float, 24>&>( st.GLoffset );
    node["GLoffset"].SetStyle(YAML::EmitterStyle::Flow);
    node["GLradius"] = st.GLradius;
    node["GLwidth"] = st.GLwidth;
    node["mode"] = st.mode;
    node["comment"] = string(st.comment);
    return node;
  };

  static bool decode(const Node& node, tpcSCGL& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    auto SC = node["SC"].as<array<float, 8>>();
    copy(begin(SC), end(SC), reinterpret_cast<float*>(st.SC));
    auto SCoffset = node["SCoffset"].as<array<float, 8>>();
    copy(begin(SCoffset), end(SCoffset), reinterpret_cast<float*>(st.SCoffset));
    auto SCexponent = node["SCexponent"].as<array<float, 8>>();
    copy(begin(SCexponent), end(SCexponent), reinterpret_cast<float*>(st.SCexponent));
    auto SCscaler = node["SCscaler"].as<array<float, 8>>();
    copy(begin(SCscaler), end(SCscaler), reinterpret_cast<float*>(st.SCscaler));
    auto GL = node["GL"].as<array<float, 24>>();
    copy(begin(GL), end(GL), reinterpret_cast<float*>(st.GL));
    auto GLoffset = node["GLoffset"].as<array<float, 24>>();
    copy(begin(GLoffset), end(GLoffset), reinterpret_cast<float*>(st.GLoffset));
    st.GLradius = node["GLradius"].as<float>();
    st.GLwidth = node["GLwidth"].as<float>();
    st.mode = node["mode"].as<int>();
    auto comment = node["comment"].as<string>();
    copy(begin(comment), end(comment), st.comment);
    st.comment[comment.size()] ='\0';
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcHVPlanes> {
  static Node encode(const tpcHVPlanes& st) {
    Node node;
    
    node["CM_shift_z"] = st.CM_shift_z;
    node["CM_tilt_x"] = st.CM_tilt_x;
    node["CM_tilt_y"] = st.CM_tilt_y;
    node["GGE_shift_z"] = st.GGE_shift_z;
    node["GGE_tilt_x"] = st.GGE_tilt_x;
    node["GGE_tilt_y"] = st.GGE_tilt_y;
    node["GGW_shift_z"] = st.GGW_shift_z;
    node["GGW_tilt_x"] = st.GGW_tilt_x;
    node["GGW_tilt_y"] = st.GGW_tilt_y;
    return node;
  };

  static bool decode(const Node& node, tpcHVPlanes& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.CM_shift_z = node["CM_shift_z"].as<float>();
    st.CM_tilt_x = node["CM_tilt_x"].as<float>();
    st.CM_tilt_y = node["CM_tilt_y"].as<float>();
    st.GGE_shift_z = node["GGE_shift_z"].as<float>();
    st.GGE_tilt_x = node["GGE_tilt_x"].as<float>();
    st.GGE_tilt_y = node["GGE_tilt_y"].as<float>();
    st.GGW_shift_z = node["GGW_shift_z"].as<float>();
    st.GGW_tilt_x = node["GGW_tilt_x"].as<float>();
    st.GGW_tilt_y = node["GGW_tilt_y"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcAltroParams> {
  static Node encode(const tpcAltroParams& st) {
    Node node;
    
    node["N"] = st.N;
    node["Altro_thr"] = st.Altro_thr;
    node["Altro_seq"] = st.Altro_seq;
    node["Altro_K1"] = st.Altro_K1;
    node["Altro_K2"] = st.Altro_K2;
    node["Altro_K3"] = st.Altro_K3;
    node["Altro_L1"] = st.Altro_L1;
    node["Altro_L2"] = st.Altro_L2;
    node["Altro_L3"] = st.Altro_L3;
    return node;
  };

  static bool decode(const Node& node, tpcAltroParams& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.N = node["N"].as<int>();
    st.Altro_thr = node["Altro_thr"].as<int>();
    st.Altro_seq = node["Altro_seq"].as<int>();
    st.Altro_K1 = node["Altro_K1"].as<int>();
    st.Altro_K2 = node["Altro_K2"].as<int>();
    st.Altro_K3 = node["Altro_K3"].as<int>();
    st.Altro_L1 = node["Altro_L1"].as<int>();
    st.Altro_L2 = node["Altro_L2"].as<int>();
    st.Altro_L3 = node["Altro_L3"].as<int>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcGridLeak> {
  static Node encode(const tpcGridLeak& st) {
    Node node;
    
    node["InnerGLRadius"] = st.InnerGLRadius;
    node["MiddlGLRadius"] = st.MiddlGLRadius;
    node["OuterGLRadius"] = st.OuterGLRadius;
    node["InnerGLWidth"] = st.InnerGLWidth;
    node["MiddlGLWidth"] = st.MiddlGLWidth;
    node["OuterGLWidth"] = st.OuterGLWidth;
    node["InnerGLStrength"] = st.InnerGLStrength;
    node["MiddlGLStrength"] = st.MiddlGLStrength;
    node["OuterGLStrength"] = st.OuterGLStrength;
    return node;
  };

  static bool decode(const Node& node, tpcGridLeak& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.InnerGLRadius = node["InnerGLRadius"].as<double>();
    st.MiddlGLRadius = node["MiddlGLRadius"].as<double>();
    st.OuterGLRadius = node["OuterGLRadius"].as<double>();
    st.InnerGLWidth = node["InnerGLWidth"].as<double>();
    st.MiddlGLWidth = node["MiddlGLWidth"].as<double>();
    st.OuterGLWidth = node["OuterGLWidth"].as<double>();
    st.InnerGLStrength = node["InnerGLStrength"].as<double>();
    st.MiddlGLStrength = node["MiddlGLStrength"].as<double>();
    st.OuterGLStrength = node["OuterGLStrength"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<g2t_vertex> {
  static Node encode(const g2t_vertex& st) {
    Node node;
    
    node["ge_volume"] = string(st.ge_volume);
    node["daughter_p"] = st.daughter_p;
    node["eg_label"] = st.eg_label;
    node["eg_proc"] = st.eg_proc;
    node["event_p"] = st.event_p;
    node["ge_medium"] = st.ge_medium;
    node["ge_proc"] = st.ge_proc;
    node["id"] = st.id;
    node["is_itrmd"] = st.is_itrmd;
    node["n_daughter"] = st.n_daughter;
    node["n_parent"] = st.n_parent;
    node["next_itrmd_p"] = st.next_itrmd_p;
    node["next_prim_v_p"] = st.next_prim_v_p;
    node["parent_p"] = st.parent_p;
    node["eg_tof"] = st.eg_tof;
    node["eg_x"] = reinterpret_cast<const array<float, 3>&>( st.eg_x );
    node["eg_x"].SetStyle(YAML::EmitterStyle::Flow);
    node["ge_tof"] = st.ge_tof;
    node["ge_x"] = reinterpret_cast<const array<float, 3>&>( st.ge_x );
    node["ge_x"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, g2t_vertex& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    auto ge_volume = node["ge_volume"].as<string>();
    copy(begin(ge_volume), end(ge_volume), st.ge_volume);
    st.ge_volume[ge_volume.size()] ='\0';
    st.daughter_p = node["daughter_p"].as<int>();
    st.eg_label = node["eg_label"].as<int>();
    st.eg_proc = node["eg_proc"].as<int>();
    st.event_p = node["event_p"].as<int>();
    st.ge_medium = node["ge_medium"].as<int>();
    st.ge_proc = node["ge_proc"].as<int>();
    st.id = node["id"].as<int>();
    st.is_itrmd = node["is_itrmd"].as<int>();
    st.n_daughter = node["n_daughter"].as<int>();
    st.n_parent = node["n_parent"].as<int>();
    st.next_itrmd_p = node["next_itrmd_p"].as<int>();
    st.next_prim_v_p = node["next_prim_v_p"].as<int>();
    st.parent_p = node["parent_p"].as<int>();
    st.eg_tof = node["eg_tof"].as<float>();
    auto eg_x = node["eg_x"].as<array<float, 3>>();
    copy(begin(eg_x), end(eg_x), reinterpret_cast<float*>(st.eg_x));
    st.ge_tof = node["ge_tof"].as<float>();
    auto ge_x = node["ge_x"].as<array<float, 3>>();
    copy(begin(ge_x), end(ge_x), reinterpret_cast<float*>(st.ge_x));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<TpcAvgCurrent> {
  static Node encode(const TpcAvgCurrent& st) {
    Node node;
    
    node["run"] = st.run;
    node["start_time"] = st.start_time;
    node["stop_time"] = st.stop_time;
    node["AvCurrent"] = reinterpret_cast<const array<float, 192>&>( st.AvCurrent );
    node["AvCurrent"].SetStyle(YAML::EmitterStyle::Flow);
    node["AcCharge"] = reinterpret_cast<const array<float, 192>&>( st.AcCharge );
    node["AcCharge"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, TpcAvgCurrent& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.run = node["run"].as<int>();
    st.start_time = node["start_time"].as<int>();
    st.stop_time = node["stop_time"].as<int>();
    auto AvCurrent = node["AvCurrent"].as<array<float, 192>>();
    copy(begin(AvCurrent), end(AvCurrent), reinterpret_cast<float*>(st.AvCurrent));
    auto AcCharge = node["AcCharge"].as<array<float, 192>>();
    copy(begin(AcCharge), end(AcCharge), reinterpret_cast<float*>(st.AcCharge));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<spaceChargeCor> {
  static Node encode(const spaceChargeCor& st) {
    Node node;
    
    node["fullFieldB"] = st.fullFieldB;
    node["halfFieldB"] = st.halfFieldB;
    node["zeroField"] = st.zeroField;
    node["halfFieldA"] = st.halfFieldA;
    node["fullFieldA"] = st.fullFieldA;
    node["satRate"] = st.satRate;
    node["factor"] = st.factor;
    node["detector"] = st.detector;
    node["offset"] = st.offset;
    node["ewratio"] = st.ewratio;
    return node;
  };

  static bool decode(const Node& node, spaceChargeCor& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.fullFieldB = node["fullFieldB"].as<double>();
    st.halfFieldB = node["halfFieldB"].as<double>();
    st.zeroField = node["zeroField"].as<double>();
    st.halfFieldA = node["halfFieldA"].as<double>();
    st.fullFieldA = node["fullFieldA"].as<double>();
    st.satRate = node["satRate"].as<double>();
    st.factor = node["factor"].as<float>();
    st.detector = node["detector"].as<float>();
    st.offset = node["offset"].as<float>();
    st.ewratio = node["ewratio"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcEffectiveGeom> {
  static Node encode(const tpcEffectiveGeom& st) {
    Node node;
    
    node["drift_length_correction"] = st.drift_length_correction;
    node["z_inner_offset"] = st.z_inner_offset;
    node["z_outer_offset"] = st.z_outer_offset;
    node["z_inner_offset_West"] = st.z_inner_offset_West;
    node["z_outer_offset_West"] = st.z_outer_offset_West;
    return node;
  };

  static bool decode(const Node& node, tpcEffectiveGeom& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.drift_length_correction = node["drift_length_correction"].as<double>();
    st.z_inner_offset = node["z_inner_offset"].as<double>();
    st.z_outer_offset = node["z_outer_offset"].as<double>();
    st.z_inner_offset_West = node["z_inner_offset_West"].as<double>();
    st.z_outer_offset_West = node["z_outer_offset_West"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<TpcResponseSimulator> {
  static Node encode(const TpcResponseSimulator& st) {
    Node node;
    
    node["I0"] = st.I0;
    node["Cluster"] = st.Cluster;
    node["W"] = st.W;
    node["OmegaTau"] = st.OmegaTau;
    node["K3IP"] = st.K3IP;
    node["K3IR"] = st.K3IR;
    node["K3OP"] = st.K3OP;
    node["K3OR"] = st.K3OR;
    node["FanoFactor"] = st.FanoFactor;
    node["AveragePedestal"] = st.AveragePedestal;
    node["AveragePedestalRMS"] = st.AveragePedestalRMS;
    node["AveragePedestalRMSX"] = st.AveragePedestalRMSX;
    node["tauIntegration"] = st.tauIntegration;
    node["tauF"] = st.tauF;
    node["tauP"] = st.tauP;
    node["tauXI"] = st.tauXI;
    node["tauXO"] = st.tauXO;
    node["tauCI"] = st.tauCI;
    node["tauCO"] = st.tauCO;
    node["SigmaJitterTI"] = st.SigmaJitterTI;
    node["SigmaJitterTO"] = st.SigmaJitterTO;
    node["SigmaJitterXI"] = st.SigmaJitterXI;
    node["SigmaJitterXO"] = st.SigmaJitterXO;
    node["longitudinalDiffusion"] = st.longitudinalDiffusion;
    node["transverseDiffusion"] = st.transverseDiffusion;
    node["NoElPerAdc"] = st.NoElPerAdc;
    node["NoElPerAdcI"] = st.NoElPerAdcI;
    node["NoElPerAdcO"] = st.NoElPerAdcO;
    node["NoElPerAdcX"] = st.NoElPerAdcX;
    node["OmegaTauScaleI"] = st.OmegaTauScaleI;
    node["OmegaTauScaleO"] = st.OmegaTauScaleO;
    node["SecRowCorIW"] = reinterpret_cast<const array<float, 2>&>( st.SecRowCorIW );
    node["SecRowCorIW"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowCorOW"] = reinterpret_cast<const array<float, 2>&>( st.SecRowCorOW );
    node["SecRowCorOW"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowCorIE"] = reinterpret_cast<const array<float, 2>&>( st.SecRowCorIE );
    node["SecRowCorIE"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowCorOE"] = reinterpret_cast<const array<float, 2>&>( st.SecRowCorOE );
    node["SecRowCorOE"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowSigIW"] = reinterpret_cast<const array<float, 2>&>( st.SecRowSigIW );
    node["SecRowSigIW"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowSigOW"] = reinterpret_cast<const array<float, 2>&>( st.SecRowSigOW );
    node["SecRowSigOW"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowSigIE"] = reinterpret_cast<const array<float, 2>&>( st.SecRowSigIE );
    node["SecRowSigIE"].SetStyle(YAML::EmitterStyle::Flow);
    node["SecRowSigOE"] = reinterpret_cast<const array<float, 2>&>( st.SecRowSigOE );
    node["SecRowSigOE"].SetStyle(YAML::EmitterStyle::Flow);
    node["PolyaInner"] = st.PolyaInner;
    node["PolyaOuter"] = st.PolyaOuter;
    node["T0offset"] = st.T0offset;
    node["T0offsetI"] = st.T0offsetI;
    node["T0offsetO"] = st.T0offsetO;
    node["FirstRowC"] = st.FirstRowC;
    return node;
  };

  static bool decode(const Node& node, TpcResponseSimulator& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.I0 = node["I0"].as<float>();
    st.Cluster = node["Cluster"].as<float>();
    st.W = node["W"].as<float>();
    st.OmegaTau = node["OmegaTau"].as<float>();
    st.K3IP = node["K3IP"].as<float>();
    st.K3IR = node["K3IR"].as<float>();
    st.K3OP = node["K3OP"].as<float>();
    st.K3OR = node["K3OR"].as<float>();
    st.FanoFactor = node["FanoFactor"].as<float>();
    st.AveragePedestal = node["AveragePedestal"].as<float>();
    st.AveragePedestalRMS = node["AveragePedestalRMS"].as<float>();
    st.AveragePedestalRMSX = node["AveragePedestalRMSX"].as<float>();
    st.tauIntegration = node["tauIntegration"].as<float>();
    st.tauF = node["tauF"].as<float>();
    st.tauP = node["tauP"].as<float>();
    st.tauXI = node["tauXI"].as<float>();
    st.tauXO = node["tauXO"].as<float>();
    st.tauCI = node["tauCI"].as<float>();
    st.tauCO = node["tauCO"].as<float>();
    st.SigmaJitterTI = node["SigmaJitterTI"].as<float>();
    st.SigmaJitterTO = node["SigmaJitterTO"].as<float>();
    st.SigmaJitterXI = node["SigmaJitterXI"].as<float>();
    st.SigmaJitterXO = node["SigmaJitterXO"].as<float>();
    st.longitudinalDiffusion = node["longitudinalDiffusion"].as<float>();
    st.transverseDiffusion = node["transverseDiffusion"].as<float>();
    st.NoElPerAdc = node["NoElPerAdc"].as<float>();
    st.NoElPerAdcI = node["NoElPerAdcI"].as<float>();
    st.NoElPerAdcO = node["NoElPerAdcO"].as<float>();
    st.NoElPerAdcX = node["NoElPerAdcX"].as<float>();
    st.OmegaTauScaleI = node["OmegaTauScaleI"].as<float>();
    st.OmegaTauScaleO = node["OmegaTauScaleO"].as<float>();
    auto SecRowCorIW = node["SecRowCorIW"].as<array<float, 2>>();
    copy(begin(SecRowCorIW), end(SecRowCorIW), reinterpret_cast<float*>(st.SecRowCorIW));
    auto SecRowCorOW = node["SecRowCorOW"].as<array<float, 2>>();
    copy(begin(SecRowCorOW), end(SecRowCorOW), reinterpret_cast<float*>(st.SecRowCorOW));
    auto SecRowCorIE = node["SecRowCorIE"].as<array<float, 2>>();
    copy(begin(SecRowCorIE), end(SecRowCorIE), reinterpret_cast<float*>(st.SecRowCorIE));
    auto SecRowCorOE = node["SecRowCorOE"].as<array<float, 2>>();
    copy(begin(SecRowCorOE), end(SecRowCorOE), reinterpret_cast<float*>(st.SecRowCorOE));
    auto SecRowSigIW = node["SecRowSigIW"].as<array<float, 2>>();
    copy(begin(SecRowSigIW), end(SecRowSigIW), reinterpret_cast<float*>(st.SecRowSigIW));
    auto SecRowSigOW = node["SecRowSigOW"].as<array<float, 2>>();
    copy(begin(SecRowSigOW), end(SecRowSigOW), reinterpret_cast<float*>(st.SecRowSigOW));
    auto SecRowSigIE = node["SecRowSigIE"].as<array<float, 2>>();
    copy(begin(SecRowSigIE), end(SecRowSigIE), reinterpret_cast<float*>(st.SecRowSigIE));
    auto SecRowSigOE = node["SecRowSigOE"].as<array<float, 2>>();
    copy(begin(SecRowSigOE), end(SecRowSigOE), reinterpret_cast<float*>(st.SecRowSigOE));
    st.PolyaInner = node["PolyaInner"].as<float>();
    st.PolyaOuter = node["PolyaOuter"].as<float>();
    st.T0offset = node["T0offset"].as<float>();
    st.T0offsetI = node["T0offsetI"].as<float>();
    st.T0offsetO = node["T0offsetO"].as<float>();
    st.FirstRowC = node["FirstRowC"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<ResponseSimulator> {
  static Node encode(const ResponseSimulator& st) {
    Node node;

    node["min_signal"] = st.min_signal;
    node["electron_cutoff_energy"] = st.electron_cutoff_energy;
    return node;
  };

  static bool decode(const Node& node, ResponseSimulator& st) {
    if(!node.IsMap()) {
      return false;
    }

    st.min_signal = node["min_signal"].as<double>();
    st.electron_cutoff_energy = node["electron_cutoff_energy"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcElectronics> {
  static Node encode(const tpcElectronics& st) {
    Node node;
    
    node["numberOfTimeBins"] = st.numberOfTimeBins;
    node["nominalGain"] = st.nominalGain;
    node["samplingFrequency"] = st.samplingFrequency;
    node["tZero"] = st.tZero;
    node["adcCharge"] = st.adcCharge;
    node["adcConversion"] = st.adcConversion;
    node["averagePedestal"] = st.averagePedestal;
    node["shapingTime"] = st.shapingTime;
    node["tau"] = st.tau;
    return node;
  };

  static bool decode(const Node& node, tpcElectronics& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.numberOfTimeBins = node["numberOfTimeBins"].as<int>();
    st.nominalGain = node["nominalGain"].as<double>();
    st.samplingFrequency = node["samplingFrequency"].as<double>();
    st.tZero = node["tZero"].as<double>();
    st.adcCharge = node["adcCharge"].as<double>();
    st.adcConversion = node["adcConversion"].as<double>();
    st.averagePedestal = node["averagePedestal"].as<double>();
    st.shapingTime = node["shapingTime"].as<double>();
    st.tau = node["tau"].as<double>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<Survey> {
  static Node encode(const Survey& st) {
    Node node;
    
    node["Id"] = st.Id;
    node["r00"] = st.r00;
    node["r01"] = st.r01;
    node["r02"] = st.r02;
    node["r10"] = st.r10;
    node["r11"] = st.r11;
    node["r12"] = st.r12;
    node["r20"] = st.r20;
    node["r21"] = st.r21;
    node["r22"] = st.r22;
    node["t0"] = st.t0;
    node["t1"] = st.t1;
    node["t2"] = st.t2;
    node["sigmaRotX"] = st.sigmaRotX;
    node["sigmaRotY"] = st.sigmaRotY;
    node["sigmaRotZ"] = st.sigmaRotZ;
    node["sigmaTrX"] = st.sigmaTrX;
    node["sigmaTrY"] = st.sigmaTrY;
    node["sigmaTrZ"] = st.sigmaTrZ;
    node["comment"] = string(st.comment);
    return node;
  };

  static bool decode(const Node& node, Survey& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.Id = node["Id"].as<int>();
    st.r00 = node["r00"].as<double>();
    st.r01 = node["r01"].as<double>();
    st.r02 = node["r02"].as<double>();
    st.r10 = node["r10"].as<double>();
    st.r11 = node["r11"].as<double>();
    st.r12 = node["r12"].as<double>();
    st.r20 = node["r20"].as<double>();
    st.r21 = node["r21"].as<double>();
    st.r22 = node["r22"].as<double>();
    st.t0 = node["t0"].as<double>();
    st.t1 = node["t1"].as<double>();
    st.t2 = node["t2"].as<double>();
    st.sigmaRotX = node["sigmaRotX"].as<double>();
    st.sigmaRotY = node["sigmaRotY"].as<double>();
    st.sigmaRotZ = node["sigmaRotZ"].as<double>();
    st.sigmaTrX = node["sigmaTrX"].as<double>();
    st.sigmaTrY = node["sigmaTrY"].as<double>();
    st.sigmaTrZ = node["sigmaTrZ"].as<double>();
    auto comment = node["comment"].as<string>();
    copy(begin(comment), end(comment), st.comment);
    st.comment[comment.size()] ='\0';
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<MagFactor> {
  static Node encode(const MagFactor& st) {
    Node node;
    
    node["ScaleFactor"] = st.ScaleFactor;
    return node;
  };

  static bool decode(const Node& node, MagFactor& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.ScaleFactor = node["ScaleFactor"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<TpcAvgPowerSupply> {
  static Node encode(const TpcAvgPowerSupply& st) {
    Node node;
    
    node["run"] = st.run;
    node["start_time"] = st.start_time;
    node["stop_time"] = st.stop_time;
    node["Current"] = reinterpret_cast<const array<float, 192>&>( st.Current );
    node["Current"].SetStyle(YAML::EmitterStyle::Flow);
    node["Charge"] = reinterpret_cast<const array<float, 192>&>( st.Charge );
    node["Charge"].SetStyle(YAML::EmitterStyle::Flow);
    node["Voltage"] = reinterpret_cast<const array<float, 192>&>( st.Voltage );
    node["Voltage"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, TpcAvgPowerSupply& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.run = node["run"].as<int>();
    st.start_time = node["start_time"].as<int>();
    st.stop_time = node["stop_time"].as<int>();
    auto Current = node["Current"].as<array<float, 192>>();
    copy(begin(Current), end(Current), reinterpret_cast<float*>(st.Current));
    auto Charge = node["Charge"].as<array<float, 192>>();
    copy(begin(Charge), end(Charge), reinterpret_cast<float*>(st.Charge));
    auto Voltage = node["Voltage"].as<array<float, 192>>();
    copy(begin(Voltage), end(Voltage), reinterpret_cast<float*>(st.Voltage));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcAnodeHV> {
  static Node encode(const tpcAnodeHV& st) {
    Node node;
    
    node["sector"] = st.sector;
    node["socket"] = st.socket;
    node["voltage"] = st.voltage;
    return node;
  };

  static bool decode(const Node& node, tpcAnodeHV& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.sector = node["sector"].as<unsigned short>();
    st.socket = node["socket"].as<unsigned short>();
    st.voltage = node["voltage"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcPadGainT0> {
  static Node encode(const tpcPadGainT0& st) {
    Node node;
    
    node["run"] = st.run;
    node["Gain"] = reinterpret_cast<const array<float, 196560>&>( st.Gain );
    node["Gain"].SetStyle(YAML::EmitterStyle::Flow);
    node["T0"] = reinterpret_cast<const array<float, 196560>&>( st.T0 );
    node["T0"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcPadGainT0& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.run = node["run"].as<int>();
    auto Gain = node["Gain"].as<array<float, 196560>>();
    copy(begin(Gain), end(Gain), reinterpret_cast<float*>(st.Gain));
    auto T0 = node["T0"].as<array<float, 196560>>();
    copy(begin(T0), end(T0), reinterpret_cast<float*>(st.T0));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<asic_thresholds> {
  static Node encode(const asic_thresholds& st) {
    Node node;
    
    node["thresh_lo"] = st.thresh_lo;
    node["thresh_hi"] = st.thresh_hi;
    node["n_seq_lo"] = st.n_seq_lo;
    node["n_seq_hi"] = st.n_seq_hi;
    return node;
  };

  static bool decode(const Node& node, asic_thresholds& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.thresh_lo = node["thresh_lo"].as<int>();
    st.thresh_hi = node["thresh_hi"].as<int>();
    st.n_seq_lo = node["n_seq_lo"].as<int>();
    st.n_seq_hi = node["n_seq_hi"].as<int>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcDimensions> {
  static Node encode(const tpcDimensions& st) {
    Node node;
    
    node["numberOfSectors"] = st.numberOfSectors;
    node["tpcInnerRadius"] = st.tpcInnerRadius;
    node["tpcOuterRadius"] = st.tpcOuterRadius;
    node["tpcTotalLength"] = st.tpcTotalLength;
    node["wheelInnerRadius"] = st.wheelInnerRadius;
    node["wheelOuterRadius"] = st.wheelOuterRadius;
    node["wheelThickness"] = st.wheelThickness;
    node["senseGasOuterRadius"] = st.senseGasOuterRadius;
    node["tpeaThickness"] = st.tpeaThickness;
    node["cathodeInnerRadius"] = st.cathodeInnerRadius;
    node["cathodeOuterRadius"] = st.cathodeOuterRadius;
    node["cathodeThickness"] = st.cathodeThickness;
    node["outerCuThickness"] = st.outerCuThickness;
    node["outerKaptonThickness"] = st.outerKaptonThickness;
    node["outerNomexThickness"] = st.outerNomexThickness;
    node["outerGlueThickness"] = st.outerGlueThickness;
    node["outerInsGasThickness"] = st.outerInsGasThickness;
    node["outerAlThickness"] = st.outerAlThickness;
    node["outerAlHoneycombThickness"] = st.outerAlHoneycombThickness;
    node["innerGlueThickness"] = st.innerGlueThickness;
    node["innerNomexThickness"] = st.innerNomexThickness;
    node["innerKaptonThickness"] = st.innerKaptonThickness;
    node["innerAlThickness"] = st.innerAlThickness;
    node["innerGapWidI"] = st.innerGapWidI;
    node["innerGapWidO"] = st.innerGapWidO;
    node["innerGapHeit"] = st.innerGapHeit;
    node["innerGapRad"] = st.innerGapRad;
    node["innerInWidth"] = st.innerInWidth;
    node["innerOutWidth"] = st.innerOutWidth;
    node["innerHeight"] = st.innerHeight;
    node["innerPPDepth"] = st.innerPPDepth;
    node["innerAlDepth"] = st.innerAlDepth;
    node["innerMWCDepth"] = st.innerMWCDepth;
    node["innerBoundary"] = st.innerBoundary;
    node["innerRCenter"] = st.innerRCenter;
    node["innerMWCInn"] = st.innerMWCInn;
    node["innerMWCOut"] = st.innerMWCOut;
    node["innerMVCHei"] = st.innerMVCHei;
    node["innerAirGaps"] = st.innerAirGaps;
    node["innerExtraAl"] = st.innerExtraAl;
    node["innerZGaps"] = reinterpret_cast<const array<double, 5>&>( st.innerZGaps );
    node["innerZGaps"].SetStyle(YAML::EmitterStyle::Flow);
    node["innerZGapsSize"] = reinterpret_cast<const array<double, 5>&>( st.innerZGapsSize );
    node["innerZGapsSize"].SetStyle(YAML::EmitterStyle::Flow);
    node["innerXExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.innerXExtraAl );
    node["innerXExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["innerZExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.innerZExtraAl );
    node["innerZExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["innerDXExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.innerDXExtraAl );
    node["innerDXExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["innerDZExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.innerDZExtraAl );
    node["innerDZExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerGapWidI"] = st.outerGapWidI;
    node["outerGapWidO"] = st.outerGapWidO;
    node["outerGapHeit"] = st.outerGapHeit;
    node["outerGapRad"] = st.outerGapRad;
    node["outerInWidth"] = st.outerInWidth;
    node["outerOutWidth"] = st.outerOutWidth;
    node["outerHeight"] = st.outerHeight;
    node["outerPPDepth"] = st.outerPPDepth;
    node["outerAlDepth"] = st.outerAlDepth;
    node["outerMWCDepth"] = st.outerMWCDepth;
    node["outerBoundary"] = st.outerBoundary;
    node["outerRCenter"] = st.outerRCenter;
    node["outerMWCInn"] = st.outerMWCInn;
    node["outerMWCOut"] = st.outerMWCOut;
    node["outerMVCHei"] = st.outerMVCHei;
    node["outerAirGaps"] = st.outerAirGaps;
    node["outerExtraAl"] = st.outerExtraAl;
    node["outerZGaps"] = reinterpret_cast<const array<double, 8>&>( st.outerZGaps );
    node["outerZGaps"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerZGapsSize"] = reinterpret_cast<const array<double, 8>&>( st.outerZGapsSize );
    node["outerZGapsSize"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerXExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.outerXExtraAl );
    node["outerXExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerZExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.outerZExtraAl );
    node["outerZExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerDXExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.outerDXExtraAl );
    node["outerDXExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    node["outerDZExtraAl"] = reinterpret_cast<const array<double, 5>&>( st.outerDZExtraAl );
    node["outerDZExtraAl"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcDimensions& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.numberOfSectors = node["numberOfSectors"].as<int>();
    st.tpcInnerRadius = node["tpcInnerRadius"].as<double>();
    st.tpcOuterRadius = node["tpcOuterRadius"].as<double>();
    st.tpcTotalLength = node["tpcTotalLength"].as<double>();
    st.wheelInnerRadius = node["wheelInnerRadius"].as<double>();
    st.wheelOuterRadius = node["wheelOuterRadius"].as<double>();
    st.wheelThickness = node["wheelThickness"].as<double>();
    st.senseGasOuterRadius = node["senseGasOuterRadius"].as<double>();
    st.tpeaThickness = node["tpeaThickness"].as<double>();
    st.cathodeInnerRadius = node["cathodeInnerRadius"].as<double>();
    st.cathodeOuterRadius = node["cathodeOuterRadius"].as<double>();
    st.cathodeThickness = node["cathodeThickness"].as<double>();
    st.outerCuThickness = node["outerCuThickness"].as<double>();
    st.outerKaptonThickness = node["outerKaptonThickness"].as<double>();
    st.outerNomexThickness = node["outerNomexThickness"].as<double>();
    st.outerGlueThickness = node["outerGlueThickness"].as<double>();
    st.outerInsGasThickness = node["outerInsGasThickness"].as<double>();
    st.outerAlThickness = node["outerAlThickness"].as<double>();
    st.outerAlHoneycombThickness = node["outerAlHoneycombThickness"].as<double>();
    st.innerGlueThickness = node["innerGlueThickness"].as<double>();
    st.innerNomexThickness = node["innerNomexThickness"].as<double>();
    st.innerKaptonThickness = node["innerKaptonThickness"].as<double>();
    st.innerAlThickness = node["innerAlThickness"].as<double>();
    st.innerGapWidI = node["innerGapWidI"].as<double>();
    st.innerGapWidO = node["innerGapWidO"].as<double>();
    st.innerGapHeit = node["innerGapHeit"].as<double>();
    st.innerGapRad = node["innerGapRad"].as<double>();
    st.innerInWidth = node["innerInWidth"].as<double>();
    st.innerOutWidth = node["innerOutWidth"].as<double>();
    st.innerHeight = node["innerHeight"].as<double>();
    st.innerPPDepth = node["innerPPDepth"].as<double>();
    st.innerAlDepth = node["innerAlDepth"].as<double>();
    st.innerMWCDepth = node["innerMWCDepth"].as<double>();
    st.innerBoundary = node["innerBoundary"].as<double>();
    st.innerRCenter = node["innerRCenter"].as<double>();
    st.innerMWCInn = node["innerMWCInn"].as<double>();
    st.innerMWCOut = node["innerMWCOut"].as<double>();
    st.innerMVCHei = node["innerMVCHei"].as<double>();
    st.innerAirGaps = node["innerAirGaps"].as<int>();
    st.innerExtraAl = node["innerExtraAl"].as<int>();
    auto innerZGaps = node["innerZGaps"].as<array<double, 5>>();
    copy(begin(innerZGaps), end(innerZGaps), reinterpret_cast<double*>(st.innerZGaps));
    auto innerZGapsSize = node["innerZGapsSize"].as<array<double, 5>>();
    copy(begin(innerZGapsSize), end(innerZGapsSize), reinterpret_cast<double*>(st.innerZGapsSize));
    auto innerXExtraAl = node["innerXExtraAl"].as<array<double, 5>>();
    copy(begin(innerXExtraAl), end(innerXExtraAl), reinterpret_cast<double*>(st.innerXExtraAl));
    auto innerZExtraAl = node["innerZExtraAl"].as<array<double, 5>>();
    copy(begin(innerZExtraAl), end(innerZExtraAl), reinterpret_cast<double*>(st.innerZExtraAl));
    auto innerDXExtraAl = node["innerDXExtraAl"].as<array<double, 5>>();
    copy(begin(innerDXExtraAl), end(innerDXExtraAl), reinterpret_cast<double*>(st.innerDXExtraAl));
    auto innerDZExtraAl = node["innerDZExtraAl"].as<array<double, 5>>();
    copy(begin(innerDZExtraAl), end(innerDZExtraAl), reinterpret_cast<double*>(st.innerDZExtraAl));
    st.outerGapWidI = node["outerGapWidI"].as<double>();
    st.outerGapWidO = node["outerGapWidO"].as<double>();
    st.outerGapHeit = node["outerGapHeit"].as<double>();
    st.outerGapRad = node["outerGapRad"].as<double>();
    st.outerInWidth = node["outerInWidth"].as<double>();
    st.outerOutWidth = node["outerOutWidth"].as<double>();
    st.outerHeight = node["outerHeight"].as<double>();
    st.outerPPDepth = node["outerPPDepth"].as<double>();
    st.outerAlDepth = node["outerAlDepth"].as<double>();
    st.outerMWCDepth = node["outerMWCDepth"].as<double>();
    st.outerBoundary = node["outerBoundary"].as<double>();
    st.outerRCenter = node["outerRCenter"].as<double>();
    st.outerMWCInn = node["outerMWCInn"].as<double>();
    st.outerMWCOut = node["outerMWCOut"].as<double>();
    st.outerMVCHei = node["outerMVCHei"].as<double>();
    st.outerAirGaps = node["outerAirGaps"].as<int>();
    st.outerExtraAl = node["outerExtraAl"].as<int>();
    auto outerZGaps = node["outerZGaps"].as<array<double, 8>>();
    copy(begin(outerZGaps), end(outerZGaps), reinterpret_cast<double*>(st.outerZGaps));
    auto outerZGapsSize = node["outerZGapsSize"].as<array<double, 8>>();
    copy(begin(outerZGapsSize), end(outerZGapsSize), reinterpret_cast<double*>(st.outerZGapsSize));
    auto outerXExtraAl = node["outerXExtraAl"].as<array<double, 5>>();
    copy(begin(outerXExtraAl), end(outerXExtraAl), reinterpret_cast<double*>(st.outerXExtraAl));
    auto outerZExtraAl = node["outerZExtraAl"].as<array<double, 5>>();
    copy(begin(outerZExtraAl), end(outerZExtraAl), reinterpret_cast<double*>(st.outerZExtraAl));
    auto outerDXExtraAl = node["outerDXExtraAl"].as<array<double, 5>>();
    copy(begin(outerDXExtraAl), end(outerDXExtraAl), reinterpret_cast<double*>(st.outerDXExtraAl));
    auto outerDZExtraAl = node["outerDZExtraAl"].as<array<double, 5>>();
    copy(begin(outerDZExtraAl), end(outerDZExtraAl), reinterpret_cast<double*>(st.outerDZExtraAl));
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcGas> {
  static Node encode(const tpcGas& st) {
    Node node;
    
    node["barometricPressure"] = st.barometricPressure;
    node["inputTPCGasPressure"] = st.inputTPCGasPressure;
    node["nitrogenPressure"] = st.nitrogenPressure;
    node["gasPressureDiff"] = st.gasPressureDiff;
    node["inputGasTemperature"] = st.inputGasTemperature;
    node["outputGasTemperature"] = st.outputGasTemperature;
    node["flowRateArgon1"] = st.flowRateArgon1;
    node["flowRateArgon2"] = st.flowRateArgon2;
    node["flowRateMethane"] = st.flowRateMethane;
    node["percentMethaneIn"] = st.percentMethaneIn;
    node["ppmOxygenIn"] = st.ppmOxygenIn;
    node["flowRateExhaust"] = st.flowRateExhaust;
    node["percentMethaneOut"] = st.percentMethaneOut;
    node["ppmWaterOut"] = st.ppmWaterOut;
    node["ppmOxygenOut"] = st.ppmOxygenOut;
    node["flowRateRecirculation"] = st.flowRateRecirculation;
    return node;
  };

  static bool decode(const Node& node, tpcGas& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.barometricPressure = node["barometricPressure"].as<float>();
    st.inputTPCGasPressure = node["inputTPCGasPressure"].as<float>();
    st.nitrogenPressure = node["nitrogenPressure"].as<float>();
    st.gasPressureDiff = node["gasPressureDiff"].as<float>();
    st.inputGasTemperature = node["inputGasTemperature"].as<float>();
    st.outputGasTemperature = node["outputGasTemperature"].as<float>();
    st.flowRateArgon1 = node["flowRateArgon1"].as<float>();
    st.flowRateArgon2 = node["flowRateArgon2"].as<float>();
    st.flowRateMethane = node["flowRateMethane"].as<float>();
    st.percentMethaneIn = node["percentMethaneIn"].as<float>();
    st.ppmOxygenIn = node["ppmOxygenIn"].as<float>();
    st.flowRateExhaust = node["flowRateExhaust"].as<float>();
    st.percentMethaneOut = node["percentMethaneOut"].as<float>();
    st.ppmWaterOut = node["ppmWaterOut"].as<float>();
    st.ppmOxygenOut = node["ppmOxygenOut"].as<float>();
    st.flowRateRecirculation = node["flowRateRecirculation"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<trgTimeOffset> {
  static Node encode(const trgTimeOffset& st) {
    Node node;
    
    node["offset"] = st.offset;
    node["laserOffset"] = st.laserOffset;
    node["laserOffsetW"] = st.laserOffsetW;
    return node;
  };

  static bool decode(const Node& node, trgTimeOffset& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    st.offset = node["offset"].as<float>();
    st.laserOffset = node["laserOffset"].as<float>();
    st.laserOffsetW = node["laserOffsetW"].as<float>();
    return true;
  }
};
}


namespace YAML {
template<>
struct convert<tpcPadrowT0> {
  static Node encode(const tpcPadrowT0& st) {
    Node node;
    
    node["T0"] = reinterpret_cast<const array<float, 100>&>( st.T0 );
    node["T0"].SetStyle(YAML::EmitterStyle::Flow);
    return node;
  };

  static bool decode(const Node& node, tpcPadrowT0& st) {
    if(!node.IsMap()) {
      return false;
    }
    
    auto T0 = node["T0"].as<array<float, 100>>();
    copy(begin(T0), end(T0), reinterpret_cast<float*>(st.T0));
    return true;
  }
};
}
