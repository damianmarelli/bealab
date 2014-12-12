/// @file bealab/core/prelim/concurrency.hpp
/// Support for concurrent programming.

#ifndef _BEALAB_PRELIM_CONCURRENCY_
#define	_BEALAB_PRELIM_CONCURRENCY_

#ifndef BEALAB_NOMPI
#include <boost/mpi.hpp>
#endif
#include <thread>
#include <mutex>
#include <future>

namespace bealab
{
/// @defgroup prelim_concurrency Concurrency
/// Support for concurrent programming.
/// @{

// Use namespace
#ifndef BEALAB_NOMPI
namespace mpi = boost::mpi;
#endif

// Use classes
using std::thread;
using std::mutex;
using std::lock_guard;
using std::unique_lock;
using std::condition_variable;

/*
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
		function<void()> task = std::bind( fun, args... );						//XXX Binding to avoid a bug in gcc when passing variadic arguments to a lambda

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
*/
#ifndef BEALAB_NOMPI

// Forward declarations
template<class value_type> class vectorx;
template<class base> class vector_interface;
template<class value_type> using Vec = vector_interface<vectorx<value_type>>;

/// Implements a for loop in parallel over a cluster.
/// It evaluates fun(i), for i=0,...,I-1, and return a vector with the results.
template<class T>
Vec<T> parallel_for( int I, const function<T(int)>& fun )
{
	// Parallel processing
	Vec<T> X(I);
	mpi::communicator com;
	#pragma omp parallel for
	for( int i = com.rank(); i < I; i += com.size() )
		X[i] = fun(i);

	// Collect the result
	vector<Vec<T>> XX;
	mpi::gather( com, X, XX, 0 );
	Vec<T> Y(I);
	if( com.rank() == 0 )
		for( int i = 0; i < I; i++ )
			Y[i] = XX[ i % com.size() ](i);

	// Share the result with the other nodes
	mpi::broadcast( com, Y, 0 );

	return Y;
}
#endif

/// @}
}
#endif
