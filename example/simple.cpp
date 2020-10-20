#include <iomanip>
#include <iostream>
#include <fstream>

#include "arg_parser.h"
#include "tpcrs/tpcrs.h"

using namespace std;


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
  simulator.Distort(begin(simu_hits), end(simu_hits), back_inserter(dist_hits));

  vector<tpcrs::DigiHit> digi_hits;
  simulator.Digitize(begin(simu_hits), end(simu_hits), back_inserter(digi_hits));

  cout << "Converted "
    << simu_hits.size() << " simulated hits into "
    << dist_hits.size() << " distorted hits and "
    << digi_hits.size() << " digitized channels\n";

  return EXIT_SUCCESS;
}
