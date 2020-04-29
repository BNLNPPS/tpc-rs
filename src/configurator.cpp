#include <iostream>
#include <sys/stat.h>

#include "tpcrs/configurator.h"
#include "config.h"

namespace tpcrs {

bool Configurator::Configure(std::string configname)
{
  name = configname;

  std::string paths(TPCRS_CONFIG_SEARCH_PATHS);

  size_t last = 0;
  size_t next = 0;
  while ((next = paths.find(':', last)) != std::string::npos)
  {
    Instance().searchPaths.push_back(paths.substr(last, next-last));
    last = next + 1;
  }
  Instance().searchPaths.push_back(paths.substr(last));

  yaml = YAML::LoadFile(Configurator::File());
}


std::string Configurator::Locate(std::string filename)
{
  struct stat buffer;   
  for (std::string path : Instance().searchPaths)
  {
    std::string fullname(path + "/" + filename);
    if (stat(fullname.c_str(), &buffer) == 0)
      return fullname;
  }

  return "";
}

}
