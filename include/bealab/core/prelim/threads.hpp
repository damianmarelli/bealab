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
