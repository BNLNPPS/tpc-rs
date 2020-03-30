#ifndef tpcrs_config_structs_h
#define tpcrs_config_structs_h

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "tpcrs/config_yaml.h"
#include "tpcrs/configurator.h"

namespace tpcrs {

struct IConfigStruct
{
  bool IsMarked() const { return isMarked; };
  void Mark(bool mark=true) { isMarked = mark; }
  void UnMark() { isMarked = false; }

  virtual unsigned int GetNRows() const { return 0; }
  virtual ~IConfigStruct() = default;
  virtual std::string GetName() const { return "noname"; }

  /// By default "Mark" this "table"/"chair" as "bad"
  bool isMarked = true;
};


template<typename Base_t, typename Chair_t, typename Struct_t>
struct ConfigStruct : Base_t
{
  static Chair_t* instance()
  {
    static Chair_t instance;
    return &instance;
  }

  Struct_t* Struct(int i=0) { return &rows_[i]; }
  Struct_t* Struct(int i=0) const { return const_cast<Struct_t*>(&rows_[i]); }
  virtual std::string GetName() const { return name; }
  unsigned int GetNRows() const { return rows_.size(); }

  static std::string name;

 protected:

  ConfigStruct()
  {
    // Deal with optionally present structs
    if (!Configurator::YAML(name)) {
      rows_.push_back(Struct_t());
      // Still create an object but "Mark" this "table"/"chair" as "bad"
      IConfigStruct::Mark();
    } else {
      try {
        rows_.push_back( Configurator::YAML(name).as< Struct_t >() );
        IConfigStruct::UnMark();
      } catch (std::exception& e) {
        rows_ = Configurator::YAML(name).as< std::vector<Struct_t> >();
        IConfigStruct::UnMark();
      }
    }
  }

  /**
   * Deleted members prohibit instantiation of this class.
   */
  ///@{
  ConfigStruct(ConfigStruct const&)   = delete;
  void operator=(ConfigStruct const&) = delete;
  ///@}

  std::vector<Struct_t> rows_;
};

}

#endif
