// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/matlab.hpp
/// Interface to call Matlab functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOMATLAB

#ifndef _BEALAB_MATLAB_
#define	_BEALAB_MATLAB_

#include <bealab/core/octave.hpp>

namespace bealab
{
/// Matlab interface
namespace matlab
{
/// @defgroup interfaces_matlab Matlab interface
/// Interface to call Matlab functions.
/// @{

// Alias for octave::functor
template<class... oargs>
using ofunctor = octave::functor<oargs...>;

/// Functor holding an Octave function
template<class... oargs>
class functor : public ofunctor<oargs...> {

	string call_prefix() override
	{
		return "matlab -nojvm -r";
	}

	string call_suffix() override
	{
		return "> /dev/null 2> /dev/null";
	}

public:
	using ofunctor<oargs...>::functor;
};

/// @}
}
}
#endif

#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
