// * $Id: StGlobalDirection.cc,v 1.1 2004/03/05 17:22:54 fisyak Exp $
#include "StDbUtilities/StGlobalDirection.hh"
static const char rcsid[] = "$Id: StGlobalDirection.cc,v 1.1 2004/03/05 17:22:54 fisyak Exp $";
std::ostream &operator<<(std::ostream &os, const StGlobalDirection &a)
{
  return os << "GC direction ( "
         << a.position().x() << ", "
         << a.position().y() << ", "
         << a.position().z() << ")";
}
