#include <iomanip>
#include <iostream>
#include <fstream>

#include "tpcrs/tpcrs.h"

using namespace std;

class ArgParser
{
 public:
  ArgParser(int &argc, char **argv);
  string get_value(const string &option) const;
  int verify() const;
 private:
   vector<string> args;
};


int main(int argc, char **argv)
{
  ArgParser arg_parser(argc, argv);

  if (arg_parser.verify())
    return EXIT_FAILURE;

  string   conf_file(arg_parser.get_value("-c"));
  ifstream data_file(arg_parser.get_value("-f"));

  tpcrs::Configurator cfg("simple", conf_file);
  tpcrs::Simulator simulator(cfg);

  vector<tpcrs::SimulatedHit> simu_hits;
  tpcrs::SimulatedHit simu_hit;

  while (data_file >> simu_hit)
    simu_hits.push_back(simu_hit);

  vector<tpcrs::DistortedHit> dist_hits;
  tpcrs::distort(begin(simu_hits), end(simu_hits), back_inserter(dist_hits), simulator);

  vector<tpcrs::DigiHit> digi_hits;
  tpcrs::digitize(begin(simu_hits), end(simu_hits), back_inserter(digi_hits), simulator);

  cout << "Converted "
    << simu_hits.size() << " simulated hits into "
    << dist_hits.size() << " distorted hits and "
    << digi_hits.size() << " digitized channels\n";

  return EXIT_SUCCESS;
}


ArgParser::ArgParser(int &argc, char **argv)
{
  for (int i=1; i < argc; ++i)
    args.push_back( string(argv[i]) );
}

string ArgParser::get_value(const string &option) const
{
  auto itr = find(args.begin(), args.end(), option);

  if (itr != args.end() && ++itr != args.end()) {
    return *itr;
  }
  return "";
}

int ArgParser::verify() const
{
  string file_name = get_value("-c");
  if (file_name.empty()) {
    cerr << "Error: Config file must be specified with -c option\n";
    return EXIT_FAILURE;
  } else if (!ifstream(file_name).is_open()) {
    cerr << "Error: Config file \"" << file_name << "\" not found\n";
    return EXIT_FAILURE;
  }

  file_name = get_value("-f");
  if (file_name.empty()) {
    cerr << "Error: Data file must be specified with -f option\n";
    return EXIT_FAILURE;
  } else if (!ifstream(file_name).is_open()) {
    cerr << "Error: Data file \"" << file_name << "\" not found\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
