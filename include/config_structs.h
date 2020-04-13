#ifndef tpcrs_config_structs_h
#define tpcrs_config_structs_h

#include <exception>
#include <string>
#include <vector>

#include "TTable.h"

#include "yaml-cpp/yaml.h"

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
};

}

#endif
