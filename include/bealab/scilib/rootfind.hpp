// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/rootfind.hpp
/// Root-finding routines.

#ifndef _BEALAB_ROOTFIND_
#define	_BEALAB_ROOTFIND_

#include <bealab/core/blas.hpp>

namespace bealab
{
/// @defgroup rootfind Root-finding
/// Root-finding routines.
/// @{

/// Finds the zero of a function within a given interval.
double fzero( function<double(double)> fun, double x_lo, double x_hi,
		int max_iter=1000, double epsabs = 1e-6, double epsrel = 1e-6 );

/// @}
}
#endif
