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
		return "matlab -nodesktop -nosplash -r";
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
