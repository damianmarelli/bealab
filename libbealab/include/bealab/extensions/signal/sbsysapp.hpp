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
