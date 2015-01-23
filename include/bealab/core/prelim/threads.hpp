// This software is licensed under the BSD 3-Clause License with the possibily to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
// You may not use this work except in compliance with the License.
// You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the file License for the specific language governing permissions and limitations under the License. 
// If you wish to obtain a commercial license, please contact the authors via e-mail.
//
// Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)

/// @file bealab/core/prelim/threads.hpp
/// Support for multithreading.

#ifndef _BEALAB_PRELIM_THREADS_
#define	_BEALAB_PRELIM_THREADS_

#include <bealab/core/prelim/imports.hpp>

namespace bealab
{
/// @defgroup prelim_threads Threads
/// Support for multithreading.
/// @{

/// Thread queue manager with a maximum number of concurrent threads.
class thread_queue {

	int max_nthreads;															///< Maximum number of concurrent threads
	int nthreads;																///< Current number of concurrent threads
	mutex qmut;																	///< Mutex for accessing the task queue and counter
	queue<function<void()>> tqueue;												///< Task queue
	condition_variable wcondvar;												///< Condition variable for waiting
	bool waiting;																///< Flag to indicate that the main thread is waiting for all threads to finish

	/// Function executed by each thread
	void thread_driver()
	{
		while(true) {

			// Retrieve a task from the queue
			function<void()> task;
			{
				lock_guard<mutex> lock(qmut);
				if( tqueue.empty() )
					break;

				task = tqueue.front();
				tqueue.pop();
			}

			// Execute the task
			task();
		}

		// On exit, decrement the thread counter and wake up the main thread if necessary
		{
			lock_guard<mutex> lock(qmut);
			nthreads--;
			if( nthreads == 0 && waiting == true )
				wcondvar.notify_all();
		}
	}

public:

	/// Constructor
	thread_queue( int N=sysconf(_SC_NPROCESSORS_ONLN) ) :
		max_nthreads(N), nthreads(0), waiting(false) {}

	/// Add a task to the queue
	template<class F, class... A>
	void async( const F& fun, const A&... args )
	{
		// Make a task
		function<void()> task = std::bind( fun, args... );

		// Enqueue the task
		{
			lock_guard<mutex> lock(qmut);
			tqueue.push( task );
		}

		// If necessary, launch a new thread
		if( nthreads < max_nthreads ) {
			lock_guard<mutex> lock(qmut);
			thread launcher( [this]() { this->thread_driver(); } );
			launcher.detach();
			nthreads++;
		}
	}

	/// Wait for all the threads to finish
	void wait()
	{
	    unique_lock<mutex> lock(qmut);
		if( nthreads == 0 )
			return;
		waiting = true;
	    wcondvar.wait( lock );
	}
};

/// @}
}
#endif
