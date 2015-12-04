/*
EAKMeans is a fast Exact K-means library written in C++ with 
command-line interface, shared library + header files and 
Python bindings

Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

This file is part of EAKMeans.

EAKMeans is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

EAKMeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.



*/

#ifndef SORTUTIL_H
#define SORTUTIL_H

#include <functional>

#include <vector>
#include <algorithm>
#include <exception>
#include <iostream>

namespace sort{

//RVO is fantastic
template <typename IndexType, typename SortableType, typename TContainer, typename SortCondition>
std::vector<IndexType> get_argsort(const TContainer & container, const SortCondition & sorter){
	std::vector<std::pair<SortableType, IndexType> > pairs;
	pairs.reserve(container.size());
	for (IndexType i = 0; i < container.size(); ++i){
		auto pair = std::make_pair(container[i], i);
		pairs.push_back(pair);
	}
	std::sort(pairs.begin(), pairs.end(), sorter);
	
	std::vector<IndexType> sorted_indices (container.size());
	for (unsigned i = 0; i < container.size(); ++i){
		sorted_indices[i] = pairs[i].second;
	}
	return sorted_indices;
}

//like above but a bit lower level
template <typename TInt, typename TSort, typename SortCondition>
void set_argsorted(TInt ndata, const TSort * const data, TInt * const argsorted, const SortCondition & sorter){
	std::vector<std::pair<TSort, TInt> > pairs;
	pairs.reserve(ndata);
	for (TInt i = 0; i < ndata; ++i){
		pairs.emplace_back( std::make_pair(data[i], i) );
	}
	std::sort(pairs.begin(), pairs.end(), sorter);
	for (TInt i = 0; i < ndata; ++i){
		argsorted[i] = pairs[i].second;
	}
}

template <typename TInt, typename TSort>
void set_argsorted_increasing(TInt ndata, const TSort * const data, TInt * const argsorted){
	set_argsorted<TInt, TSort, std::less<std::pair<TSort, TInt> > > (ndata, data, argsorted, std::less<std::pair<TSort, TInt> >());
}

template <typename IndexType, typename SortableType, typename TContainer>
std::vector<IndexType> get_argsort_decreasing(const TContainer & container){
	return get_argsort<IndexType, SortableType, TContainer, std::greater<std::pair<SortableType, IndexType> > > (container, std::greater<std::pair<SortableType, IndexType> >());
}


template <typename IndexType, typename SortableType, typename TContainer>
std::vector<IndexType> get_argsort_increasing(const TContainer & container){
	return get_argsort<IndexType, SortableType, TContainer, std::less<std::pair<SortableType, IndexType> > > (container, std::less<std::pair<SortableType, IndexType> >());
}




/* Sort a Container (an array for example) into left and right partitions, according to the bool flags (true) left (false) right provided in isleft_flag. Input -- n_left: how many flags say true. stride: how many objects (of type *TIterator) does a flag correspond to. swap_function: how should two memory sections of objects be swapped. */
template <typename TIterator, typename TSwapFunction>
void partition_left_right(const TIterator & a, unsigned n_left, unsigned stride, const bool * const isleft_flags, const TSwapFunction & swap_function){
	TIterator left = a; 
	const TIterator b = a + stride*n_left;
	TIterator right = a + stride*n_left;
	const bool * left_flag = isleft_flags;
	const bool * right_flag = isleft_flags + n_left;
	
	/* while left partition has not finished being scanned for roght objects, */
	while (left != b){
		
		/* find next bad left */
		while (left !=  b && *left_flag == true){
			++left_flag;
			left += stride;
		}
	
		if (left == b){
			return;
		}
	
		/* find next bad right */
		while (*right_flag == false){
			++right_flag;
			right += stride;
		}
		
		/* swap bad left with bad right */
		swap_function(stride, left, right);
		right += stride;
		left += stride;
		++left_flag;
		++right_flag;
	}
}




/* for N items, std::nth_element returns partitions an array into items less than nth largest and items greater than nth largest. This function does the same kind of thing, but for several ns, where the ns are 2**k for k = 0:(less than npartitions). Due to the application at hand, this function foes geometric_elements for nrows of size ncols 
 * Examples, if npartitions is 1
 * x x x x x x x x x x x x x x x x x x x x 
 * if npartitions is 2:
 * x o x x x x x x x x x x x x x x x x x x 
 * if npartitions is 3: 
 * x o x o x x x x x x x x x x x x x x x x 
 * if npartitions is 4:
 * x o x o x x x o x x x x x x x x x x x x
 * if npartitions is 5:
 * x o x o x x x o x x x x x x x o x x x x  
 * 
 * where the o are in the correct position, as if completely sorted.
 *  */
template <typename IndexType, typename SortableType, typename SortCondition>
void geometric_elements(IndexType nrows, IndexType ncols, IndexType npartitions, SortableType * const data, const SortCondition & sorter){
	
	if (ncols <= std::pow(2, npartitions - 1)){
		throw std::runtime_error("too many partitions requested in sorting algorithm set_geometric_elements. 2^(npartitions - 1) cannot exceed ncols");
	}
	
	
	
	IndexType middleindex;
	IndexType endindex;
	
	for (IndexType r = 0; r < nrows;  ++r){
		middleindex = ncols;
		for (IndexType pi = npartitions - 1; pi > 0; --pi){
			endindex = middleindex;
			middleindex = std::pow(2, pi) - 1;
			std::nth_element(data + r*ncols, data + r*ncols + middleindex, data + r*ncols + endindex, sorter); 
		} 
	}
}

/* 1-D version of geomteric_elements (without test of ncols vs npartitions) */
template <typename IndexType, typename SortableType, typename SortCondition>
void geometric_elements(IndexType ncols, IndexType npartitions, SortableType * const data, const SortCondition & sorter){
	IndexType middleindex = ncols;
	IndexType endindex;	
	for (IndexType pi = npartitions - 1; pi > 0; --pi){
		endindex = middleindex;
		middleindex = std::pow(2, pi) - 1;
		std::nth_element(data, data + middleindex, data + endindex, sorter); 
	} 
}



		

//used in annulus-like kmeans. I call it exponion k-means
template <typename TInt, typename TSortable, typename SortCondition>
void update_halfordered_partvals_etc1(TInt ncentroids, TInt nrows, TInt nparts, const TSortable * const CC, std::pair<TSortable, TInt> * const halforderedpairs, TSortable * const partvals, TInt * const indices_halforderedCC, const SortCondition & sorter){
	if (nparts != static_cast<TInt>(std::floor(std::log2(static_cast<double>(ncentroids - 1))))){
		throw std::runtime_error("I was of the opinion that npartitions would be constrained to be std::floor(std::log2(static_cast<double>(ncentroids - 1)) in this function, has something gone awry?");
	}

	//we do not force nrows to be ncentroids for mulithreading purposes
	TInt size_indices_halforderedCC = 2*(std::pow(2, nparts - 1) - 1);
	TInt delta_part;
	TInt sortindices_start;
	TInt final_relevant_index = 2;
	for (TInt i = 0; i<nparts -2; ++i){
		final_relevant_index*=2;
	}
	
	for (TInt r = 0; r < nrows; ++r){
		for (TInt c= 0; c < ncentroids; ++ c){
			halforderedpairs[r*ncentroids + c].first = CC[r*ncentroids + halforderedpairs[r*ncentroids + c].second]*0.5;
		}
		//TODO: use geometric_elements 1-d version
		std::partial_sort(halforderedpairs + r*ncentroids, halforderedpairs + r*ncentroids + final_relevant_index, halforderedpairs + (r + 1)*ncentroids, sorter);
		delta_part = 1;
		sortindices_start = 0;
		for (TInt i = 0; i < nparts - 1; ++i){
			delta_part *= 2;
			partvals[r*(nparts - 1) + i] = halforderedpairs[ncentroids*r + delta_part - 1].first;
			for (TInt j = 0; j<delta_part; ++j){
				indices_halforderedCC[r*size_indices_halforderedCC  + sortindices_start + j] = halforderedpairs[r*ncentroids + j].second;
			}
			std::sort(indices_halforderedCC + r*size_indices_halforderedCC  + sortindices_start, indices_halforderedCC + r*size_indices_halforderedCC  + sortindices_start + delta_part);
			
			sortindices_start += delta_part;
		}
	}	
} 


//used in r12v3
template <typename TInt, typename TSortable, typename SortCondition>
void update_halfordered_partvals_nosort1(TInt ncentroids, TInt nrows, TInt nparts, const TSortable * const CC, std::pair<TSortable, TInt> * const halforderedpairs, TSortable * const partvals, TInt * const indices_halforderedCC, const SortCondition & sorter){
	if (nparts != static_cast<TInt>(std::floor(std::log2(static_cast<double>(ncentroids - 1))))){
		throw std::runtime_error("I was of the opinion that npartitions would be constrained to be std::floor(std::log2(static_cast<double>(ncentroids - 1)) in this function, has something gone awry?");
	}

	//we do not force nrows to be ncentroids for mulithreading purposes
	TInt size_indices_halforderedCC = 2*(std::pow(2, nparts - 1) - 1);
	TInt delta_part;
	TInt sortindices_start;
	TInt final_relevant_index = 2;
	for (TInt i = 0; i<nparts -2; ++i){
		final_relevant_index*=2;
	}
	
	for (TInt r = 0; r < nrows; ++r){
		for (TInt c= 0; c < ncentroids; ++ c){
			halforderedpairs[r*ncentroids + c].first = CC[r*ncentroids + halforderedpairs[r*ncentroids + c].second]*0.5;
		}
		//TODO: use geometric_elements 1-d version
		std::partial_sort(halforderedpairs + r*ncentroids, halforderedpairs + r*ncentroids + final_relevant_index, halforderedpairs + (r + 1)*ncentroids, sorter);
		delta_part = 1;
		sortindices_start = 0;
		for (TInt i = 0; i < nparts - 1; ++i){
			delta_part *= 2;
			partvals[r*(nparts - 1) + i] = halforderedpairs[ncentroids*r + delta_part - 1].first;
			for (TInt j = 0; j<delta_part; ++j){
				indices_halforderedCC[r*size_indices_halforderedCC  + sortindices_start + j] = halforderedpairs[r*ncentroids + j].second;
			}
			
			//only difference from etc1
			//std::sort(indices_halforderedCC + r*size_indices_halforderedCC  + sortindices_start, indices_halforderedCC + r*size_indices_halforderedCC  + sortindices_start + delta_part);
			
			sortindices_start += delta_part;
		}
	}	
}


//used in r12v4
template <typename TInt, typename TSortable, typename SortCondition>
void update_halfordered_partvals_nosort2(TInt ncentroids, TInt nrows, TInt nparts, const TSortable * const CC, std::pair<TSortable, TInt> * const halforderedpairs, TSortable * const partvals, TInt * const indices_halforderedCC, const SortCondition & sorter){
	if (nparts != static_cast<TInt>(std::floor(std::log2(static_cast<double>(ncentroids - 1))))){
		throw std::runtime_error("I was of the opinion that npartitions would be constrained to be std::floor(std::log2(static_cast<double>(ncentroids - 1)) in this function, has something gone awry?");
	}

	//we do not force nrows to be ncentroids for mulithreading purposes
	TInt size_indices_halforderedCC = 2*(std::pow(2, nparts - 1) - 1);
	TInt delta_part;
	TInt sortindices_start;
	TInt final_relevant_index = 2;
	for (TInt i = 0; i<nparts -2; ++i){
		final_relevant_index*=2;
	}
	
	for (TInt r = 0; r < nrows; ++r){
		for (TInt c= 0; c < ncentroids; ++ c){
			halforderedpairs[r*ncentroids + c].first = CC[r*ncentroids + halforderedpairs[r*ncentroids + c].second]*0.5;
		}


		geometric_elements(static_cast<TInt>(1), ncentroids, nparts, halforderedpairs + r*ncentroids, sorter);
		delta_part = 1;
		sortindices_start = 0;
		for (TInt i = 0; i < nparts - 1; ++i){
			delta_part *= 2;
			partvals[r*(nparts - 1) + i] = halforderedpairs[ncentroids*r + delta_part - 1].first;
			for (TInt j = 0; j<delta_part; ++j){
				indices_halforderedCC[r*size_indices_halforderedCC  + sortindices_start + j] = halforderedpairs[r*ncentroids + j].second;
			}
			
			//only difference from etc1
			//std::sort(indices_halforderedCC + r*size_indices_halforderedCC  + sortindices_start, indices_halforderedCC + r*size_indices_halforderedCC  + sortindices_start + delta_part);
			
			sortindices_start += delta_part;
		}
	}	
}



//see comments at top of r12v5 for clarification as to what this function does
template <typename TInt, typename TSortable, typename SortCondition>
void update_pairs_parts_indices_halvies(TInt nrows, TInt ncols, TInt npartitions, const TSortable * const CC, 
std::pair<TSortable, TInt> * const geometricpairs_halvies, TSortable * const partitionvalues_halvies, TInt * const geometricindices, const SortCondition & sorter){
	
	TInt delta_part;
	
	TInt ncm1 = ncols - 1;
	for (TInt r = 0; r < nrows; ++r){
		//reset the float values in the pairs from CC
		for (TInt c= 0; c < ncm1; ++ c){
			geometricpairs_halvies[r*ncm1 + c].first = CC[r*ncols + geometricpairs_halvies[r*ncm1 + c].second]*0.5;
		}
		/* perform geometric sorting on the pairs so that they are
		 * (nearest, [2nd], 3rd, [4th, 5th, 6th], 7th, [8th, 9th, ... 14th], 15th, etc)
		 * */
		geometric_elements(ncm1, npartitions, geometricpairs_halvies + r*ncm1, sorter);
		
		/* populate partvals: 
		 * extract the pivot values from geometricpairs_halvies, the pivot values being at indices [0,2,6,14...]
		 * */
		delta_part = 1;
		for (TInt i = 0; i < npartitions - 1; ++i){
			delta_part *= 2;
			partitionvalues_halvies[r*(npartitions - 1) + i] = geometricpairs_halvies[ncm1*r + delta_part - 2].first;
		}
		
		/* populate geometricindices: 
		 * */
		for (TInt c= 0; c < ncm1; ++ c){
			geometricindices[r*ncm1 + c] = geometricpairs_halvies[r*ncm1 + c].second;
		}
	
		//if (r == 5){
			//std::cout << "\n\ngeometricpairs_halvies : " << std::endl;
			//for (TInt c= 0; c < ncm1; ++ c){
				//std::cout << "(" << geometricpairs_halvies[r*ncm1 + c].first << " " << geometricpairs_halvies[r*ncm1 + c].second << ") ";
			//}
			//std::cout << std::endl;
			//std::abort();
		//}
	}	
}





}

#endif
