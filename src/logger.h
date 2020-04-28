#ifndef TPCRS_LOGGER_H_
#define TPCRS_LOGGER_H_

#include <iostream>

extern std::ostream& LOG_INFO;
extern std::ostream& LOG_WARN;
extern std::ostream& LOG_ERROR;
extern std::ostream& LOG_FATAL;
extern std::ostream& LOG_DEBUG;

std::ostream& endm(std::ostream& os);

#endif
