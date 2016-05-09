/// @file bealab/extensions/other/mlearn.hpp
/// Machine learning methods.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOGCLASSES

#ifndef _BEALAB_MLEARN_
#define _BEALAB_MLEARN_

#include <bealab/scilib.hpp>

namespace bealab
{
/// Machine learning module
namespace machine_learning
{
/// @defgroup mlearn Machine learning
/// Machine learning methods.
/// @{

class neural_network {

	void* pnn;
	void* prng;
	int dim_in  = 0;
	int dim_out = 0;

public:

	neural_network( const ivec& layer_nodes );

	~neural_network();

	void train( const vec<rvec>& input, const vec<rvec>& output );

	rvec predict( const rvec& in );

	void save( const string& fname, const string& vname="w" );

	void load( const string& fname, const string& vname="w" );
};

/// @}
}
}
#endif
#endif
