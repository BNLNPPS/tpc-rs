#pragma once

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "tpcrs/configurator.h"
#include "tpcrs/config_yaml.h"


namespace tpcrs {

struct IConfigStruct
{
  void Initialize() {}
  virtual unsigned int GetNRows() const { return 0; }
  virtual ~IConfigStruct() = default;
  virtual std::string GetName() const { return "noname"; }

  /// By default mark this structure as missing
  bool is_missing = true;
};


template<typename Base_t, typename Chair_t, typename Struct_t>
struct ConfigStruct : Base_t
{
  Struct_t* Struct(int i=0) { return &rows_[i]; }
  Struct_t* Struct(int i=0) const { return const_cast<Struct_t*>(&rows_[i]); }
  virtual std::string GetName() const { return name; }
  unsigned int GetNRows() const { return rows_.size(); }

  static std::string name;

  friend Configurator;

 protected:

  ConfigStruct(const Configurator& cfg) : rows_(), cfg_(cfg)
  {
    // Deal with optionally present structs
    if (!cfg_.YAML(name)) {
      rows_.push_back(Struct_t());
      Base_t::is_missing = true;
    } else {
      try {
        rows_.push_back( cfg_.YAML(name).as< Struct_t >() );
      } catch (std::exception& e) {
        rows_ = cfg_.YAML(name).as< std::vector<Struct_t> >();
      }
      Base_t::is_missing = false;
    }

    Base_t::Initialize();
  }

  void operator=(const ConfigStruct &other)
  {
    rows_ = other.rows_;
    const_cast<Configurator&>(cfg_) = other.cfg_;
  }

  std::vector<Struct_t> rows_;

  const Configurator& cfg_;
};

}
