/// @file bealab/extensions/signal/sparse.hpp
/// Base class for sparse approximation algorithms.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SPARSE_BASE_
#define	_BEALAB_EXTENSIONS_SIGNAL_SPARSE_BASE_

#include <bealab/core/blas.hpp>
#include <bealab/core/plot.hpp>
#include <bealab/scilib/linalg.hpp>
#include <bealab/scilib/optimization.hpp>
#include <bealab/scilib/calculus.hpp>

namespace bealab
{
namespace signal
{
namespace sparse
{
/// @defgroup sparse-base Base
/// Base class for sparse approximation algorithms.
/// @{

/// Base for all sparse approximation algorithms (abstract class).
/// @param RAN range type
/// @param COE coefficient type
/// @param IDX index type
template<class COE=double>
class base {
protected:

	/// Function for plotting
	virtual
	void plotfun() {}

public:

	/// Virtual destructor
	virtual
	~base() {}

	/// @name Parameters
	bool plot        = false;													///< Plot flag
	bool trace       = false;													///< Tracing flag
	int max_sparsity = imax;													///< Maximum number of non-zero entries (-1 = unbounded)
	double tolerance = 0;														///< Approximation error tolerance
	/// @}

    /// @name Solution
    ivec indexes;																///< Indexes of the sparse solution vector
	Vec<COE> coefficients;														///< Values (at the indexes) of the sparse solution vector
	double error;																///< Final approximation error
	/// @}

	/// Do the sparse approximation
	virtual
	void approximate() = 0;
};

/// @}
}
}
}
#endif
