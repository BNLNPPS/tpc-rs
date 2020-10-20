#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>


class ArgParser
{
 public:
  ArgParser(int &argc, char **argv);
  std::string get_value(const std::string &option) const;
  int verify();
 private:
  std::vector<std::string> args;
};


ArgParser::ArgParser(int &argc, char **argv)
{
  std::string prg_name(argv[0]);
  prg_name.erase(0, prg_name.find_last_of('/')+1);

  args.push_back(prg_name);

  for (int i=1; i < argc; ++i)
    args.push_back( std::string(argv[i]) );
}

std::string ArgParser::get_value(const std::string &option) const
{
  auto itr = std::find(args.begin(), args.end(), option);

  if (itr != args.end() && std::next(itr) != args.end() && (*++itr)[0] != '-') {
    return *itr;
  }
  return "";
}

int ArgParser::verify()
{
  std::string file_name = get_value("-c");

  if (file_name.empty()) {
    file_name = args[0] + ".yaml";
    std::cerr << "Warning: Using default config file \"" << file_name << "\". Use -c option to override\n";

    args.erase( std::remove(begin(args), end(args), "-c"), args.end() );
    args.push_back("-c");
    args.push_back(file_name);
  }

  if (!std::ifstream(file_name).is_open()) {
    std::cerr << "Error: Config file \"" << file_name << "\" not found\n";
    return EXIT_FAILURE;
  }

  file_name = get_value("-f");

  if (file_name.empty()) {
    file_name = args[0] + ".dat";
    std::cerr << "Warning: Using default data file \"" << file_name << "\". Use -f option to override\n";

    args.erase( std::remove(begin(args), end(args), "-f"), args.end());
    args.push_back("-f");
    args.push_back(file_name);
  }

  if (!std::ifstream(file_name).is_open()) {
    std::cerr << "Error: Data file \"" << file_name << "\" not found\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
