/// @file bealab/core/prelim/imports.hpp
/// Preliminary file.

#ifndef _BEALAB_PRELIM_IMPORTS_
#define	_BEALAB_PRELIM_IMPORTS_

#include <initializer_list>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <queue>
#include <set>
#include <map>
#include <algorithm>
#include <bitset>
#include <tuple>
#include <functional>
#include <exception>

namespace bealab
{
/// @defgroup prelim_imports Imports
/// Import some symbols from the STD library:
/// - Type traints
///   - declval
///   - result_of
///   - enable_if
///   - conditional
///   - is_base_of
///   - is_same
///   - is_convertible
///   - is_fundamental
/// - Input / output
///   - cin
///   - cout
///   - cerr
///   - endl
/// - Streams
///   - ostringstream
///   - ofstream
/// - Containers
///   - string
///   - vector
///   - deque
///   - queue
///   - set
///   - map
/// - Algorithms
///   - sort
/// - Other
///   - tuple
///   - get
///   - initializer_list
///   - bitset
///   - function
///   - exception
/// - Threads
///   - thread
///   - mutex
///   - lock_guard
///   - unique_lock
///   - condition_variable
/// @{

// Type traits
using std::declval;
using std::result_of;
using std::enable_if;
using std::conditional;
using std::is_base_of;
using std::is_same;
using std::is_convertible;
using std::is_fundamental;

// Input/output
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

// Streams
using std::istringstream;
using std::ostringstream;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::fstream;

// Containers
using std::string;
using std::vector;
using std::deque;
using std::queue;
using std::set;
using std::map;

// Algorithms
using std::sort;

// Other
using std::tuple;
using std::get;
using std::initializer_list;
using std::bitset;
using std::function;
using std::exception;

/// @}
}
#endif
