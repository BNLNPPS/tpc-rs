/*******************************************************************
 * Author: brian Feb 6, 1998
 * Description:  Raw data information along with access functions
 *******************************************************************/
#include "coords/StTpcPadCoordinate.hh"

// Non-Member function
std::ostream &operator<<(std::ostream &os, const StTpcPadCoordinate &a)
{
  return os << "(sector= " << a.sector()
         << ", row= "    << a.row()
         << ", pad= "    << a.pad()
         << ", tbuck= "  << a.timeBucket() << ")";
}
