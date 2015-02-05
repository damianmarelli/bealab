/// @file bealab/core/matlab.hpp
/// Interface to call Matlab functions.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOMATLAB

#ifndef _BEALAB_MATLAB_
#define	_BEALAB_MATLAB_

#include <bealab/core/octave.hpp>

namespace bealab
{
namespace matlab
{

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
		return "' > /dev/null 2> /dev/null";
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
