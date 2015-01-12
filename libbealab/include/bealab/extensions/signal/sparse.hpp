/// @file bealab/extensions/signal/sparse.hpp
/// Algorithms for sparse approximations.

#ifndef _BEALAB_EXTENSIONS_SPARSE_
#define	_BEALAB_EXTENSIONS_SPARSE_

/// @defgroup sparse Sparse approximation algorithms
/// Algorithms for sparse approximations.
/// @{

namespace bealab
{
namespace signal
{
/// Sparse approximation module
namespace sparse {}
}
}

/// @defgroup sparse-base
#include <bealab/extensions/signal/sparse/base.hpp>

/// @defgroup sparse-greedy
#include <bealab/extensions/signal/sparse/greedy.hpp>

/// @defgroup sparse-relaxation
#include <bealab/extensions/signal/sparse/relaxation.hpp>

/// @}
#endif
