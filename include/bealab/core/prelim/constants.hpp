// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
