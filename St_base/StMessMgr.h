#ifndef StMessMgr_h
#define StMessMgr_h

#include <iostream>

extern std::ostream& LOG_INFO;
extern std::ostream& LOG_WARN;
extern std::ostream& LOG_ERROR;
extern std::ostream& LOG_FATAL;
extern std::ostream& LOG_DEBUG;

std::ostream& endm(std::ostream& os);

#endif

// $Id: StMessMgr.h,v 1.13 2009/06/22 22:36:02 fine Exp $
