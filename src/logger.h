#ifndef tpcrs_logging_h
#define tpcrs_logging_h

#include <iostream>

extern std::ostream& LOG_INFO;
extern std::ostream& LOG_WARN;
extern std::ostream& LOG_ERROR;
extern std::ostream& LOG_FATAL;
extern std::ostream& LOG_DEBUG;

std::ostream& endm(std::ostream& os);

#endif
