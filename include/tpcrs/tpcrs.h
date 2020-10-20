#pragma once

#include "tpcrs/tpcrs_core.h"
#include "tpcrs/detail/digitizer.h"
#include "tpcrs/detail/distorter.h"
#include "tpcrs/detail/mag_field.h"
#include "tpcrs/detail/simulator.h"

namespace tpcrs {


class MagField : detail::MagField
{
 public:

  MagField(const tpcrs::Configurator& cfg) : detail::MagField(cfg) {}

  template<typename Vec3>
  Vec3 ValueAt(Vec3 p) const { return detail::MagField::ValueAt(p); }
};


class Distorter : detail::Distorter
{
 public:

  Distorter(const tpcrs::Configurator& cfg) : detail::Distorter(cfg) {}

  template<typename Vec3, typename MagField>
  Vec3 Distort(Vec3 p, int sector, const MagField& mag_field) const
  {
    return detail::Distorter::Distort(p, sector, mag_field);
  }
};


class Simulator : detail::Simulator
{
 public:

  Simulator(const tpcrs::Configurator& cfg) : detail::Simulator(cfg) {}

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized) const
  {
    return detail::Simulator::Digitize(first_hit, last_hit, digitized);
  }

  template<typename InputIt, typename OutputIt, typename MagField>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized, const MagField& mag_field) const
  {
    return detail::Simulator::Digitize(first_hit, last_hit, digitized, mag_field);
  }

  template<typename InputIt, typename OutputIt>
  OutputIt Distort(InputIt first_hit, InputIt last_hit, OutputIt distorted) const
  {
    return detail::Simulator::Distort(first_hit, last_hit, distorted);
  }

  template<typename InputIt, typename OutputIt>
  OutputIt Simulate(InputIt first_hit, InputIt last_hit, OutputIt charges) const
  {
    return detail::Simulator::Simulate(first_hit, last_hit, charges);
  }
};


class Digitizer : detail::Digitizer
{
 public:

  Digitizer(const tpcrs::Configurator& cfg) : detail::Digitizer(cfg) {}

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(unsigned sector, InputIt first_channel, InputIt last_channel, OutputIt digitized) const
  {
    return detail::Digitizer::Digitize(sector, first_channel, last_channel, digitized);
  }
};

}
