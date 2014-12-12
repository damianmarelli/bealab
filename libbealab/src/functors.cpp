#include <bealab/scilib/functors.hpp>
#include <bealab/scilib/rootfind.hpp>

namespace bealab
{

function<double(double)> inv( const function<double(double)>& fun, double lo, double hi )
{
	return [=]( double x )
	{
		return fzero( [&fun,x](double y){ return fun(y)-x; }, lo, hi );
	};
}

}
