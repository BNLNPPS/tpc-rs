#include "coords/StGlobalDirection.hh"
std::ostream &operator<<(std::ostream &os, const StGlobalDirection &a)
{
  return os << "GC direction ( "
         << a.position().x() << ", "
         << a.position().y() << ", "
         << a.position().z() << ")";
}
