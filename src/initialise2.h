/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef INITIALISE2_H
#define INITIALISE2_H

//#include "randomarray.h"
#include <random>
#include "sample.h"
#include "arrutilv2l0.h"
#include "arrutilv2l2.h"


namespace kmeans{
namespace initialise2{

//get copyindices guaranteeing that all distinct.
template <typename TFloat, typename TInt>
std::tuple<std::unique_ptr<TFloat []>, std::vector<TInt> > get_initialisation_indices(TInt ncentroids, TInt ndata, TInt dimension, const TFloat * const data){
	
	std::vector<TInt> initialisation_indices (ncentroids);
	std::unique_ptr<TFloat []> C(new TFloat [ncentroids*dimension]);
	
	TInt nattempts = 0;
	TInt currentindex = 0;
	
	while (currentindex < ncentroids && nattempts < 5*ncentroids){
		TInt proposal = rand() % ndata;
		bool rejected = false;
		for (TInt ci = 0; ci < currentindex; ++ci){
			TFloat l1normdiff = 0;
			for (TInt d = 0; d < dimension; ++d){
				l1normdiff += std::abs(data[proposal*dimension + d] - C[ci*dimension + d]);
			}
			if (l1normdiff < 1e-5){
				//std::cout << "---+--- NTY ---+---" << std::flush;
				rejected = true;
				break;
			}
		}
		++nattempts;
		if (rejected == false){
			std::memcpy(C.get() + currentindex*dimension, data + proposal*dimension, sizeof(TFloat)*dimension);
			initialisation_indices[currentindex] = proposal;
			++currentindex;
		}
	}
	
	if (currentindex != ncentroids){
		throw std::runtime_error("Tried to find a set of distinct datapoints, but failed (nattempts/ncentroids = 5)"); 
	}
	

	return std::make_tuple (std::move(C), std::move(initialisation_indices));
}

template <typename TFloat, typename TInt> //class URNG ?
std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, TFloat > get_kmeanspp_initialisation(TInt ndata, TInt dimension, const TFloat * const data, const TFloat * const data_l22s, TInt ncentroids){
	
	
	std::cout << "in get_kmeanspp_initialisation with (ndata, dimension) : ( " << ndata << ", " << dimension << ")" << std::endl; 
	
	TFloat l22;
	TInt index;
	TFloat cum_val;
	
	std::unique_ptr<TFloat []> C_uptr (new TFloat [ncentroids*dimension]);
	std::unique_ptr<TFloat []> C_l22s_uptr (new TFloat [ncentroids]);
	std::unique_ptr<TInt []> ind0 (new TInt [ncentroids]);
	
	auto C = C_uptr.get();
	auto C_l22s = C_l22s_uptr.get();
	
	
	std::unique_ptr<TFloat []> min_r2 (new TFloat [ndata]);
	std::unique_ptr<TFloat []> cum_min_r2 (new TFloat [ndata]);
	
	for (TInt i = 0; i < ndata; ++i){
		min_r2[i] = std::numeric_limits<TFloat>::max();
	}
	
	//first centroid chosen at random.
	index = rand()%ndata;
	ind0[0] = index;
	std::memcpy(C, data + index*dimension, sizeof(TFloat)*dimension);
	arrutilv2::set_l22(dimension, C, C_l22s[0]);
	
	//all subsequent centroids chosen with p \propto r^2  where r is min(...)
	for (TInt ci = 0; ci < ncentroids - 1; ++ci){
		//datapoint index 0:
		arrutilv2::set_l22(dimension, data, C + ci*dimension, data_l22s[0], C_l22s[ci], l22); 
		if (l22 < 0){
			std::string errm = "The value of l22 (distance squared from data[0] to centroid ci) is negative, bailing";	
			throw std::logic_error(errm);			
		}
		
		min_r2[0] = std::min(min_r2[0], l22);
		cum_min_r2[0] = min_r2[0];
		//all other datapoints,
		for (TInt i = 1; i < ndata; ++i){
			arrutilv2::set_l22(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], l22); 
			min_r2[i] = std::min(min_r2[i], l22);
			cum_min_r2[i] = cum_min_r2[i - 1] +  min_r2[i];
			
			if (l22 < 0){
				std::string errm = "The value of l22 (distance squared from data[i] to centroid ci) is negative, bailing";	
				errm = errm + "\ni : " + std::to_string(i) + " ci : " + std::to_string(ci) + " data_l22s[i] : " + std::to_string(data_l22s[i]) + " C_l22s[ci] : " + std::to_string(C_l22s[ci]) +"\n";
				throw std::logic_error(errm);			
			}

		}
	
		cum_val = cum_min_r2[ndata - 1] * ( rand() / static_cast<TFloat>(RAND_MAX) );
		//determine insertion index:
		index = arrutilv2::get_insertion_index_binarial(cum_val, ndata, cum_min_r2.get());
		if (index == ndata){
			std::string errm = "The index returned is ndata : the value to insert is larger than the total cumulative value. ";
			errm = errm + "cum_val = " + std::to_string(cum_val) + " and cum_min_r2[ndata - 1] = " + std::to_string(cum_min_r2[ndata - 1]);	
			throw std::logic_error(errm);
		}
		
		ind0[ci +1] = index;
		std::memcpy(C + (ci+1)*dimension, data + index*dimension, sizeof(TFloat)*dimension);
		arrutilv2::set_l22(dimension, C + (ci+1)*dimension, C_l22s[ci+1]);
		
	}
	
	TFloat mse2 = 0;
	TInt ci = ncentroids - 1;
	for (TInt i = 0; i < ndata; ++i){
		arrutilv2::set_l22(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], l22); 
		min_r2[i] = std::min(min_r2[i], l22);
		mse2 += min_r2[i];
	}
	mse2 /= static_cast<TFloat> (ndata);
		
	return std::make_tuple< 
	std::unique_ptr<TFloat [] >, 
	std::unique_ptr<TFloat [] >,
	std::unique_ptr<TInt [] >,
	TFloat > 
	(std::move(C_uptr), std::move(C_l22s_uptr), std::move(ind0), TFloat(mse2));	
}


}
}
#endif

