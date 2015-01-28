// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
	vec<COE> coefficients;														///< Values (at the indexes) of the sparse solution vector
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
