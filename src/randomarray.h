/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef RANDOMARRAY_H
#define RANDOMARRAY_H


#include <exception>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <random>
#include <memory>

namespace randomutil{ 
namespace randomarray{
	
template <typename IntType, typename SizeType>
void filluniform_int(SizeType size_tofill, IntType * tofill, IntType lower, IntType upper){
	IntType range = upper - lower;
	for (SizeType i = 0; i < size_tofill; ++ i){
		tofill[i] = lower + rand() % range;
	}
}

template <typename FloatType, typename SizeType>
void filluniform_float(SizeType size_tofill, FloatType * tofill, FloatType lower, FloatType upper){
	FloatType range = upper - lower;
	for (SizeType i = 0; i < size_tofill; ++ i){
		tofill[i] = lower + range * (static_cast <FloatType> (rand()) / static_cast <FloatType> (RAND_MAX));
	}
}

template <typename NumberType, typename SizeType>
void filluniform(SizeType size_tofill, NumberType * tofill, NumberType lower, NumberType upper);

template<typename SizeType>
void filluniform(SizeType size_tofill, float * tofill, float lower, float upper){
	filluniform_float(size_tofill, tofill, lower, upper);
}

template<typename SizeType>
void filluniform(SizeType size_tofill, double * tofill, double lower, double upper){
	filluniform_float(size_tofill, tofill, lower, upper);
}


template<typename SizeType>
void filluniform(SizeType size_tofill, unsigned * tofill, unsigned lower, unsigned upper){
	filluniform_int(size_tofill, tofill, lower, upper);
}

template<typename SizeType>
void filluniform(SizeType size_tofill, int * tofill, int lower, int upper){
	filluniform_int(size_tofill, tofill, lower, upper);
}

/* fill tofill with values chosen from options uniformly at random */
template <typename NumberType, typename SizeType, typename Container>
void filluniform(SizeType size_tofill, NumberType * tofill, Container && options){
	unsigned n_options = options.size();
	std::vector<decltype(n_options)> option_numbers (size_tofill);
	filluniform(size_tofill, option_numbers.data(), static_cast<unsigned> (0), n_options);
	for (SizeType i = 0; i < size_tofill; ++i){
		tofill[i] = options[option_numbers[i]];
	}
}


template <typename NumberType, typename SizeType, typename Container>
std::vector<NumberType> getuniform(SizeType N, Container && options){
	std::vector<NumberType>  sampled (N,0);
	filluniform(N, sampled.data(), options);
	return sampled;
}


// untested function:
template <typename NumberType, typename SizeType>
std::vector<NumberType> getuniform(SizeType N, NumberType lower, NumberType upper){
	std::vector<NumberType>  sampled (N,0);
	filluniform(N, sampled.data(), lower, upper);
	return sampled;
}

// untested function:
template <typename NumberType, typename SizeType>
std::unique_ptr<NumberType [] > getuniform_uptr(SizeType N, NumberType lower, NumberType upper){
	std::unique_ptr<NumberType [] > sampled (new NumberType [N]);
	filluniform(N, sampled.get(), lower, upper);
	return sampled;
}


//return vector of length ndraes of sorted vectors, each vector has probability that a value (TInt) lies in the vector being p in range [0, N) and 0 otherwise.
template <typename TInt, typename TFloat>
std::vector<std::vector<TInt>> get_p_sample(TInt ndraws, TFloat p, TInt N){
	
	TInt proposal;
	bool goodproposal;
	std::vector<std::vector<TInt>> samples (ndraws);
	
	std::default_random_engine generator(rand());
	std::binomial_distribution<TInt> distribution(N,p);
	for (TInt draw=0; draw<ndraws; ++draw) {
		//number of distinct values in [0, N) to insert into samples[draw]:
    TInt number = distribution(generator);
    samples[draw].resize(number);
    TInt i = 0;    
    while (i < number){
			proposal = rand()%N;
			goodproposal = true;
			//confirm that proposal is good:
			for (TInt j = 0; j < i; ++j){
				if (proposal == samples[draw][j]){
					goodproposal = false;
					break;
				}
			}
			if (goodproposal == true){
				samples[draw][i] = proposal;
				++i;
			}
		}
		std::sort(samples[draw].begin(), samples[draw].end());
  }
  
  return samples;
}



}
}

namespace randomutil{ 
namespace noise{

template <typename NumberType, typename SizeType, typename ProbType>
void signflip(SizeType ndata, NumberType * const data, ProbType switch_probability){
	for (SizeType i = 0; i < ndata; ++i){
		if ((float(rand()) / float(RAND_MAX)) < switch_probability){
			data[i]*=(-1);
		}
	}
}
	
}
}	

#endif
