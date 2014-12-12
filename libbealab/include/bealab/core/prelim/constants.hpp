/// @file bealab/core/prelim/constants.hpp
/// Some numeric constants

#ifndef _BEALAB_PRELIM_CONSTANTS_
#define	_BEALAB_PRELIM_CONSTANTS_

#include <limits>
#include <bealab/core/prelim/complex.hpp>

namespace bealab
{
/// @defgroup prelim_constants Constants
/// Some numeric constants
/// @{

const complex i       = complex(0,1);
const complex j       = complex(0,1);
const double pi       = 3.14159265358979323846;
const double e        = 2.71828182845904523536;
const double eps      = std::numeric_limits<double>::epsilon();
const double nan      = std::numeric_limits<double>::quiet_NaN();
const double inf      = std::numeric_limits<double>::infinity();
const int imax        = std::numeric_limits<int>::max();
const char* const tab = "\t";

/// @}
}
#endif
