/// @file bealab/extensions/other/machine_learning.hpp
/// Machine learning methods.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOGCLASSES

#ifndef _BEALAB_MACHINE_LEARNING_
#define _BEALAB_MACHINE_LEARNING_

#include <bealab/scilib.hpp>

namespace bealab
{
namespace machine_learning
{
/// @defgroup machine_learning Machine learning
/// Machine learning methods.
/// @{

class neural_network {

	void* pnn;
	void* prng;
	int dim_in  = 0;
	int dim_out = 0;

public:

	neural_network( const ivec& layer_nodes );

	void train( const Vec<rvec>& input, const Vec<rvec>& output );

	rvec predict( const rvec& in );
};

/// @}
}
}
#endif
#endif
