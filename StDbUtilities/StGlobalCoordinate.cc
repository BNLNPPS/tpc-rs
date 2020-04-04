/*************************************************************************
 * Author:  brian May 20, 1998
 *
 * Description:  Raw data information along with access functions
 *************************************************************************/

#include "StDbUtilities/StGlobalCoordinate.hh"

// Non-member functions
std::ostream &operator<<(std::ostream &os, const StGlobalCoordinate &a)
{
  return os << "GC ( "
         << a.position().x() << ", "
         << a.position().y() << ", "
         << a.position().z() << ")";
}
