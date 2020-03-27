#ifndef tpcrs_Configurator_h
#define tpcrs_Configurator_h

#include <string>
#include <vector>


namespace tpcrs {

/**
 * A singleton class 
 */
class Configurator
{
 public:

  static Configurator& Instance()
  {
    static Configurator instance;
    return instance;
  }

  bool Configure(std::string cfgname = "");

  static std::string Locate(std::string filename = "");
  static std::string File() { return Locate(Instance().name + ".yaml"); };

 private:

  /**
   * Private deleted constructors prohibit any instantiation of this class.
   */
  ///@{ 
  Configurator() {}
  Configurator(Configurator const&)   = delete;
  void operator=(Configurator const&) = delete;
  ///@}

  /// A unique name used in various file names associated with this Configurator.
  std::string name;

  std::vector<std::string> searchPaths;
};

}

#endif
