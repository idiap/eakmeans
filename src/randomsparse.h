/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef RANDOMSPARSE_H
#define RANDOMSPARSE_H


#include <exception>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <random>
#include <memory>
#include "sparsedatasets.h"
#include "randomarray.h"

namespace randomutil{ 
namespace randomsparse{

//TODO : make function accept a random number generator. 
template <typename TInt, typename TFloat>
sparse::SparseData<TInt, TFloat> get_sparsedata(TInt ndata, TInt dimension, TFloat sparsity){
	
	
	std::vector<TFloat> values;
	std::vector<TInt> indices;
	std::vector<TInt> starts (1,0);
	std::vector<std::string> labels;
	
	
	//vector of vectors of indices
	auto vvinds = randomutil::randomarray::get_p_sample<TInt, TFloat>(ndata, sparsity, dimension);
	
	
	for (TInt i = 0; i < ndata; ++i){
		for (auto & index : vvinds[i]){
			indices.push_back(index);
			TFloat value = (static_cast <TFloat> (rand()) / static_cast <TFloat> (RAND_MAX));
			values.push_back(value);
		}
		//for (TInt j = 0; j < dimension; ++j){ //really slow way to do it! TODO: 
			//bool nonzero = (sparsity > (rand() / (RAND_MAX + 0.)));
			//if (nonzero == true){
				//TFloat value = (static_cast <TFloat> (rand()) / static_cast <TFloat> (RAND_MAX));
				//values.push_back(value);
				//indices.push_back(j);
			//}
		//}
		starts.push_back(indices.size());
		labels.push_back("1011"); //give everything label 0
	}
	
	return sparse::SparseData<TInt, TFloat>(std::move(values), std::move(indices), std::move(starts), std::move(labels));
	
}

//TODO : make function accept a random number generator. 
void write_sparsedata(unsigned ndata, unsigned dimension, double sparsity, const std::string & filename, bool dimheader = true){
	auto sd  = get_sparsedata<unsigned, double> (ndata, dimension, sparsity);
	sd.write(filename, dimheader); 
}

void write_sparse_and_dense_data(unsigned ndata, unsigned dimension, double sparsity, const std::string & sparsefilename,  const std::string & densefilename){
	auto sd  = get_sparsedata<unsigned, double> (ndata, dimension, sparsity);
	sd.write(sparsefilename, true); 
	sd.write_dense(densefilename, true); 
}


}

}


#endif

