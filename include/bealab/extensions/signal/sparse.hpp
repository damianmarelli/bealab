// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
