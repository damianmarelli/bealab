// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/scilib/optimization.hpp
/// Numeric parameter optimization routines.

#ifndef _BEALAB_OPTIMIZATION_
#define	_BEALAB_OPTIMIZATION_

#include <bealab/core/blas.hpp>

/// @defgroup optimization Optimization
/// Parameter optimization routines.
/// @{

/// @defgroup optimization-base
#include <bealab/scilib/optimization/base.hpp>

/// @defgroup optimization-fromscratch
#include <bealab/scilib/optimization/fromscratch.hpp>

/// @defgroup optimization-gsl
#include <bealab/scilib/optimization/gsl.hpp>

/// @defgroup optimization-nlopt
#include <bealab/scilib/optimization/nlopt.hpp>

/// @defgroup optimization-ipopt
#include <bealab/scilib/optimization/ipopt.hpp>

/// @}
#endif
