// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/extensions/other/mlearn.hpp
/// Machine learning methods.

#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOGCLASSES

#ifndef _BEALAB_MLEARN_
#define _BEALAB_MLEARN_

#include <bealab/scilib.hpp>

namespace bealab
{
/// Machine learning module
namespace machine_learning
{
/// @defgroup mlearn Machine learning
/// Machine learning methods.
/// @{

class neural_network {

	void* pnn;
	void* prng;
	int dim_in  = 0;
	int dim_out = 0;

public:

	neural_network( const ivec& layer_nodes );

	void train( const Vec<rvec>& input, const Vec<rvec>& output );

	rvec predict( const rvec& in );
};

/// @}
}
}
#endif
#endif
