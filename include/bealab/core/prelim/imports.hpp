// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

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
#include <thread>
#include <mutex>
#include <future>

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
/// - Threads
///   - thread
///   - mutex
///   - lock_guard
///   - unique_lock
///   - condition_variable
/// - Tuples
///   - tuple
///   - get
///   - tie
///   - tuple_cat
/// - Other
///   - initializer_list
///   - bitset
///   - function
///   - exception
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

// Threads
using std::thread;
using std::mutex;
using std::lock_guard;
using std::unique_lock;
using std::condition_variable;

// Tuples
using std::tuple;
using std::get;
using std::tie;
using std::tuple_cat;

// Other
using std::initializer_list;
using std::bitset;
using std::function;
using std::exception;

/// @}
}
#endif
