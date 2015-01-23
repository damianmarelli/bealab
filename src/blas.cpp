// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

#include <bealab/core/blas.hpp>

namespace bealab
{

rvec vslice( double from, double skip, double N )
{
//	int N = abs( floor( (to-from) / skip ) + 1 );
	rvec x(N);
	for( int n = 0; n < N; n++ )
		x(n) = from + n*skip;
	return x;
}

ivec vrange( int from, int to_plus_1 )
{
	return vslice( from, 1, to_plus_1-from );
}

rvec linspace( double from, double to, int N )
{
//	if( N < 0 ) {
//		assert(  mod(to-from,1) == 0 );
//		N = abs(to-from) + 1;
//	}

	// Case N == 0
	if( N == 0 )
		return rvec();

	// Case from == to
	if( from == to )
		return from * ones(N);

	// Case N == 1 and from ~= to
	if( N == 1 )
		error("linspace() - At least two points need to be generated");

	// Normal case: N > 1 and from ~= to
	double skip = (to-from) / (N-1);
	rvec x(N);
	for( int n = 0; n < N; n++ )
		x(n) = from + n*skip;
	return x;
}

rvec logspace( double a, double b, int N )
{
	double alog = log(a);
	double blog = log(b);
	rvec xlog   = linspace( alog, blog, N );
	return exp(xlog);
}

}
