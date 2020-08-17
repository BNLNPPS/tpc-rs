#pragma once

#include <bitset>

namespace tpcrs { namespace detail {


template<typename Enum, bool IsEnum = std::is_enum<Enum>::value>
class BitSet;


template<typename Enum>
class BitSet<Enum, true>
{
 public:

  BitSet() = default;
  BitSet(Enum value) : bits(1 << static_cast<std::size_t>(value)) {}
  BitSet(const BitSet& other) : bits(other.bits) {}

  BitSet operator| (Enum value) const { BitSet result = *this; result.bits |= 1 << static_cast<std::size_t>(value); return result; }
  BitSet operator& (Enum value) const { BitSet result = *this; result.bits &= 1 << static_cast<std::size_t>(value); return result; }
  BitSet operator^ (Enum value) const { BitSet result = *this; result.bits ^= 1 << static_cast<std::size_t>(value); return result; }
  BitSet operator~ () const { BitSet result = *this; result.bits.flip(); return result; }

  BitSet& operator|= (Enum value) { bits |= 1 << static_cast<std::size_t>(value); return *this; }
  BitSet& operator&= (Enum value) { bits &= 1 << static_cast<std::size_t>(value); return *this; }
  BitSet& operator^= (Enum value) { bits ^= 1 << static_cast<std::size_t>(value); return *this; }

  bool any() const { return bits.any(); }
  bool all() const { return bits.all(); }
  bool none() const { return bits.none(); }
  operator bool() const { return any(); }

  bool test(Enum value) const { return bits.test(1 << static_cast<std::size_t>(value)); }
  void set(Enum value) { bits.set(1 << static_cast<std::size_t>(value)); }
  void unset(Enum value) { bits.reset(1 << static_cast<std::size_t>(value)); }

 private:

  std::bitset<64> bits;
};


template<typename Enum>
typename std::enable_if<std::is_enum<Enum>::value, BitSet<Enum>>::type
operator| (Enum left, Enum right)
{
  return BitSet<Enum>(left) | right;
}


template<typename Enum>
typename std::enable_if<std::is_enum<Enum>::value, BitSet<Enum>>::type
operator& (Enum left, Enum right)
{
  return BitSet<Enum>(left) & right;
}


template<typename Enum>
typename std::enable_if<std::is_enum<Enum>::value, BitSet<Enum>>::type
operator^ (Enum left, Enum right)
{
  return BitSet<Enum>(left) ^ right;
}


} }
