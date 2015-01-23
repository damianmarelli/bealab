// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/signal/sbsysapp.hpp
/// Subband method for system approximation.

#ifndef _BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_
#define	_BEALAB_EXTENSIONS_SIGNAL_SBSYSAPP_

/// @defgroup sbsysapp Subband system approximation
/// Subband method for system approximation.
/// @{

namespace bealab
{
namespace signal
{
/// Subband system approximation module
namespace sbsysapp {}
}
}

/// @defgroup sbsysapp_sbconfig
#include <bealab/extensions/signal/sbsysapp/sbconfig.hpp>

/// @defgroup sbsysapp_sbapprox
#include <bealab/extensions/signal/sbsysapp/sbapprox.hpp>

/// @defgroup sbsysapp_sbgreedy
#include <bealab/extensions/signal/sbsysapp/sbgreedy.hpp>

/// @defgroup sbsysapp_sbrelax
#include <bealab/extensions/signal/sbsysapp/sbrelax.hpp>

/// @}
#endif
