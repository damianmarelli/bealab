/*******************************************************************************
 * This software is licensed under the BSD 3-Clause License with the possibility to obtain a commercial license, if you cannot abide by the terms of the BSD 3-Clause license.
 * You may not use this work except in compliance with the License.
 * You may obtain a copy of the License at: http://opensource.org/licenses/BSD-3-Clause
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the file License for the specific language governing permissions and limitations under the License. 
 * If you wish to obtain a commercial license, please contact the authors via e-mail.
 *
 * Copyright (c) 2015, Damian Marelli (Damian.Marelli@newcastle.edu.au)
 *******************************************************************************/
#include <bealab/core/prelim/config.hpp>
#ifndef BEALAB_NOGCLASSES

#include <bealab/extensions/other/mlearn.hpp>
#include <GClasses/GNeuralNet.h>
#include <GClasses/GActivation.h>
#include <GClasses/GRand.h>
#include <GClasses/GDom.h>
#include <GClasses/GLearner.h>

using namespace GClasses;

namespace bealab
{
namespace machine_learning
{
//------------------------------------------------------------------------------
// Class neural_network
//------------------------------------------------------------------------------
neural_network::neural_network( const ivec& layer_nodes )
{
	// Init
	prng             = new GRand(0);
	pnn              = new GNeuralNet(*reinterpret_cast<GRand*>(prng));
	GNeuralNet* pnnr = reinterpret_cast<GNeuralNet*>(pnn);

	// Set hidden layers
	pnnr->setActivationFunction(new GActivationLogistic(), true);
//	pnnr->setActivationFunction(new GActivationArcTan(), true);
//	pnnr->setActivationFunction(new GActivationGaussian(), true);
	int L = layer_nodes.size();
	for( int l = 0; l < L; l++ )
		pnnr->addLayer( layer_nodes(l) );

	// Set output layer
//	pnnr->setActivationFunction(new GActivationIdentity(), true);
	pnnr->setActivationFunction(new GActivationBiDir(), true);
}

neural_network::~neural_network()
{
	delete reinterpret_cast<GRand*>(prng);
	delete reinterpret_cast<GNeuralNet*>(pnn);
}

void neural_network::train( const vec<rvec>& input, const vec<rvec>& output )
{
	assert( input.size() == output.size() );

	// Input data
	int I  = input.size();
	dim_in = input(0).size();
	GMatrix in(I,dim_in);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < dim_in; j++ )
			in[i][j] = input(i)(j);

	// output data
	dim_out = output(0).size();
	GMatrix out(I,dim_out);
	for( int i = 0; i < I; i++ )
		for( int j = 0; j < dim_out; j++ )
			out[i][j] = output(i)(j);

	// Train the NN
	reinterpret_cast<GNeuralNet*>(pnn)->train( in, out );
}

rvec neural_network::predict( const rvec& in )
{
	rvec out(dim_out);
	reinterpret_cast<GNeuralNet*>(pnn)->predict( &in(0), &out(0) );
	return out;
}

void neural_network::save( const string& fname, const string& vname )
{
//	int N = reinterpret_cast<GNeuralNet*>(pnn)->countWeights();
//	rvec w(N);
//	reinterpret_cast<GNeuralNet*>(pnn)->weights( &w(0) );
//	matfile mf(fname);
//	mf.save( vname, w );
//	mf.save( "dim_in", dim_in );
//	mf.save( "dim_out", dim_out );
	GDom doc;
	doc.setRoot( reinterpret_cast<GNeuralNet*>(pnn)->serialize(&doc) );
	doc.saveJson("my_neural_net.json");

}

void neural_network::load( const string& fname, const string& vname )
{
//	// Load data
//	matfile mf(fname);
//	dim_in  = mf.load<int>( "dim_in" );
//	dim_out = mf.load<int>( "dim_out" );
//	rvec w  = mf.load<rvec>( vname );
//
//	// Initialize NN
////	GMatrix inputs( 0, dim_in );
////	GMatrix outputs( 0, dim_out );
////	reinterpret_cast<GNeuralNet*>(pnn)->beginIncrementalLearning( inputs.relation(), outputs.relation() );
//	reinterpret_cast<GNeuralNet*>(pnn)->setWeights( &w(0) );

//	GRand p(0);
//	GDom doc;
//	doc.loadJson("my_neural_net.json");
//	GLearnerLoader ll(p);
//	reinterpret_cast<GNeuralNet*>(pnn) = (GNeuralNet*)ll.loadLearner(doc.root());
}

}
}
#endif

