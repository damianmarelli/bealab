/// @file bealab/extensions/other/networks.hpp
/// Generic classes for building networks (e.g., for consensus).

#ifndef _BEALAB_NETWORKS_
#define _BEALAB_NETWORKS_

#include <bealab/core/blas.hpp>

namespace bealab
{
namespace networks
{
/// @defgroup networks Networks
/// Generic classes for building networks (e.g., for consensus).
/// @{

// Forward declaration
template<class message> class node_base;

/// Link class.
/// Models a link with template messages.
template<class message>
class link {

	message inbuff;																///< Link's input
	message outbuff;															///< Link's output

public:

	typedef message message_t;													///< Message type

	node_base<message>* source;													///< Pointer to the link's source node
	node_base<message>* range;													///< Pointer to the link's range node

	/// Put a message in the link
	void write( const message& val ) { inbuff = val; }

	/// Get a message from the link
	message read() { return outbuff; }

	/// Process the link (i.e., pass the input message to the output)
	void trigger() { outbuff = inbuff; };
};

/// Base class for all node types.
/// Models a node with template messages.
template<class message>
class node_base {
public:

	typedef message message_t;													///< Message type

	vector<link<message>*> inputs;												///< Input links of the node
	vector<link<message>*> outputs;												///< Output links of the node

	/// Virtual destructor
	virtual
	~node_base() {}

	///< Process one step
	virtual
	void trigger() = 0;
};

/// Generic directed graph.
/// Models an array of template nodes, connected with directed links.
template<class node>
class directed_graph {
public:

	typedef typename node::message_t message_t;									///< Node's message type

	vector<node> nodes;															///< Nodes of the network
	vector<link<message_t>> links;												///< Links of the network

	/// Constructor
	directed_graph( const deque<node>& initialized_nodes, const deque<ivec>& source_range_map ) :
		nodes( initialized_nodes.begin(), initialized_nodes.end() )
	{
		// Connect nodes with links
		int L = source_range_map.size();
		links.resize(L);
		for( int l = 0; l < L; l++ ) {
			int source = source_range_map[l](0);
			int range  = source_range_map[l](1);
			nodes[source].outputs.push_back( &links[l] );
			nodes[range ].inputs. push_back( &links[l] );
			links[l].source = &nodes[source];
			links[l].range  = &nodes[range];
		}
	}

	/// Make the network evolve one step
	void trigger()
	{
		// Process nodes
		int N = nodes.size();
		for( int n = 0; n < N; n++ )
			nodes[n].trigger();

//		// Process nodes (with threads)
//		thread_queue tq;
//		int N = nodes.size();
//		for( int n = 0; n < N; n++ )
//			tq.async( [n,this]() { this->nodes[n].trigger(); } );
//		tq.wait();


		// Process links
		int L = links.size();
		for( int l = 0; l < L; l++ )
			links[l].trigger();

//		// Process links (with threads)
//		int L = links.size();
//		for( int l = 0; l < L; l++ )
//			tq.async( [l,this]() { this->links[l].trigger(); } );
//		tq.wait();
	}
};

/// @}
}
}
#endif
