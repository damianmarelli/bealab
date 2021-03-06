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

using bealab::octave::_function;

/// Functor handler for a Matlab function
template<class... oargs>
class function : public _function<oargs...> {

	string prefix() override
	{
		return "matlab -nojvm -r '";
	}

	/// Suffix to the shell command
	string suffix() override
	{
		return "' > /dev/null";
	}

public:

	using _function<oargs...>::_function;
};

/// @}
}
}
#endif

#endif

/*
 * TI VOGLIO TANTO BENE :* (smak!)
 */
