// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

#include <bealab/core/blas/prelim.hpp>

namespace bealab
{

void error( const string &str )
{
	cerr << str << endl;
	abort();
}

void warning( const string &str )
{
	cerr << str << endl;
}

int system( const string& call )
{
	int rv = ::system( call.data() );
//	if( rv != 0 )
//		warning("bealab::system() - Error calling: " + call);
	return rv;
}

int rw_pipe( const string& program, const vector<string>& arguments,
		int* write_fd, int* read_fd, int* error_fd )
{
	// Init
	const int READ   = 0;
	const int WRITE  = 1;
	const int STDIN  = 0;
	const int STDOUT = 1;
	const int STDERR = 2;

	// Parent -> child pipe
	int pc_pipe[2];
	if( pipe( pc_pipe ) )
		return -1;

	// Child -> parent pipe
	int cp_pipe[2];
	if( pipe( cp_pipe ) )
		return -1;

	// Child -> parent pipe (error channel)
	int cp_err_pipe[2];
	if( pipe( cp_err_pipe ) )
		return -1;

	// Calling parameters
	int P = arguments.size();
	const char* args[P+2];
	for( int p = 0; p < P; p++ )
		args[p+1] = arguments[p].data();
	args[0] = program.data();
	args[P+1] = NULL;

	// Fork
	pid_t pid = fork();
	if( pid == -1 )
		return -1;

	// Child process
	if( !pid ) {

		// Close unused pipe ends
		close( pc_pipe[WRITE] );
		close( cp_pipe[READ] );
		close( cp_err_pipe[READ] );

		// Redirect standard I/O to pipes
		dup2( pc_pipe[READ], STDIN );
		dup2( cp_pipe[WRITE], STDOUT );
		dup2( cp_err_pipe[WRITE], STDERR );

		// Start maxima
		execv( program.data(), const_cast<char* const*>(args) );
		return -1;
	}
	// Parent process
	else {

		// Close unused pipe ends
		close( pc_pipe[READ] );
		close( cp_pipe[WRITE] );
		close( cp_err_pipe[WRITE] );

		// Take the IDs of the pipe
		*write_fd = pc_pipe[WRITE];
		*read_fd  = cp_pipe[READ];
		if( error_fd )
			*error_fd = cp_err_pipe[READ];
	}

	// Return the result
	return 0;
}
}
