// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal.hpp
/// Includes all headers from signal processing modules

#ifndef _BEALAB_EXTENSIONS_SIGNAL_HPP_
#define	_BEALAB_EXTENSIONS_SIGNAL_HPP_

namespace bealab
{
/// Signal processing module
namespace signal {}
}

/// @defgroup signal Signal processing
/// @{

/// @defgroup mrsigpro
#include <bealab/extensions/signal/mrsigpro.hpp>

/// @defgroup tfanalysis
#include <bealab/extensions/signal/tfanalysis.hpp>

/// @defgroup sparse
#include <bealab/extensions/signal/sparse.hpp>

/// @defgroup sbsysapp
#include <bealab/extensions/signal/sbsysapp.hpp>

/// @defgroup misc
#include <bealab/extensions/signal/misc.hpp>

/// @}
#endif
