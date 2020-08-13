#pragma once

#include "tpcrs/tpcrs_core.h"
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


class Simulator : detail::Simulator
{
 public:

  Simulator(const tpcrs::Configurator& cfg) : detail::Simulator(cfg) {}

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(InputIt first_hit, InputIt last_hit, OutputIt digitized)
  {
    return detail::Simulator::Digitize(first_hit, last_hit, digitized);
  }

  template<typename InputIt, typename OutputIt>
  OutputIt Distort(InputIt first_hit, InputIt last_hit, OutputIt distorted)
  {
    return detail::Simulator::Distort(first_hit, last_hit, distorted);
  }
};

}
