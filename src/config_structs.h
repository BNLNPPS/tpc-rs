#ifndef TPCRS_CONFIG_STRUCTS_H_
#define TPCRS_CONFIG_STRUCTS_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "tpcrs/configurator.h"
#include "config_yaml.h"


namespace tpcrs {

struct IConfigStruct
{
  bool IsMarked() const { return isMarked; };
  void Mark(bool mark=true) { isMarked = mark; }
  void UnMark() { isMarked = false; }

  void Initialize() {}
  virtual unsigned int GetNRows() const { return 0; }
  virtual ~IConfigStruct() = default;
  virtual std::string GetName() const { return "noname"; }

  /// By default "Mark" this "table"/"chair" as "bad"
  bool isMarked = true;
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

  ConfigStruct(const Configurator& cfg) : cfg_(cfg)
  {
    // Deal with optionally present structs
    if (!cfg_.YAML(name)) {
      rows_.push_back(Struct_t());
      // Still create an object but "Mark" this "table"/"chair" as "bad"
      IConfigStruct::Mark();
    } else {
      try {
        rows_.push_back( cfg_.YAML(name).as< Struct_t >() );
        IConfigStruct::UnMark();
      } catch (std::exception& e) {
        rows_ = cfg_.YAML(name).as< std::vector<Struct_t> >();
        IConfigStruct::UnMark();
      }
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

#endif
