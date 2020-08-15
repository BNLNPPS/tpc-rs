#pragma once

#include "tpcrs/configurator.h"
#include "tpcrs/detail/enum_bitset.h"


namespace tpcrs { namespace detail {


class Distorter
{
 public:

  enum class Distortions : unsigned {
    kNone = 0,
    kMagneticField
  };

  Distorter(const Configurator& cfg, Distortions distortions = Distortions::kNone) :
    cfg_(cfg),
    distortions_(distortions)
  {}

  template<typename Vec3>
  Vec3 Distort(Vec3 p) const;

 private:

  const tpcrs::Configurator& cfg_;

  BitSet<Distortions> distortions_;

  /** Distortion due to the shape of the magnetic field */
  template<typename Vec3>
  Vec3 DistortDueMagneticField(Vec3 p) const;
};



template<typename Vec3>
Vec3 Distorter::Distort(Vec3 p) const
{
  Vec3 p_prime{};

  if ( distortions_.test(Distortions::kMagneticField) ) {
    p_prime = DistortDueMagneticField(p);
  }

  return p_prime;
}


template<typename Vec3>
Vec3 Distorter::DistortDueMagneticField(Vec3 p) const
{
  return p;
}

} }
