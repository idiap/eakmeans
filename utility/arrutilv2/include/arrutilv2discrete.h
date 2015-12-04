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

#ifndef ARRUTILV2DISCRETE_H
#define ARRUTILV2DISCRETE_H

#include <memory>
#include <vector>
#include <limits>

/* histograms, where non zero, where above threshold etc. 
 * templates should :
 * -- have no typename TFloat (better off in arrutilv2lxxx)*/

//TODO : think about when to return std::unique_ptr< []> and when to return std::vector

namespace arrutilv2{
 
template <typename TInt, typename TLabel>
std::unique_ptr<TLabel []> get_vhisto(TInt nrows, TInt ncols, TInt range, const TLabel * const labels){
	std::unique_ptr<TInt []> vhisto (new TLabel [range*ncols]);
	std::fill_n(vhisto.get(), range*ncols, 0);
	for (TInt r = 0; r < nrows; ++r){
		for (TInt c = 0; c < ncols; ++ c){
			++vhisto[ncols*labels[r*ncols + c] + c];
		}
	}
	return vhisto; 
}

template <typename TInt, typename TLabel>
std::vector<TInt> get_where_nonzero(TInt N, const TLabel * const X){
	std::vector<TInt> where;
	for (TInt i = 0; i < N; ++i){
		if (X[i] > 0){
			where.push_back(i);
		}
	}
	return where;
}


template <typename TInt, typename TLabel>
std::vector<TInt> get_where_above_threshold(TInt N, const TLabel * const X, TInt threshold){
	std::vector<TInt> where;
	for (TInt i = 0; i < N; ++i){
		if (X[i] > threshold){
			where.push_back(i);
		}
	}
	return where;
}




template <typename TInt>
std::unique_ptr<TInt []> gethistogram(TInt N, TInt range, const TInt * const L){
	std::unique_ptr<TInt []> hist(new TInt[range]); 
	std::fill_n(hist.get(), range, 0);
	for (TInt i = 0; i < N; ++i){
		++hist[L[i]];
	}
	return hist;
}

/* set max in partitions where partitions are ~equal  */
template <typename TInt, typename TNumber>
void set_maxinpartition(TInt npts, TInt npartitions, const TNumber * const vals, TNumber * const maxinpartition){
	TInt p_end = 0;
	TInt p_start;
	for (TInt p = 0; p < npartitions; ++p){
		p_start = p_end;
		p_end = ((p+1)*npts)/npartitions;
		maxinpartition[p] = std::numeric_limits<TNumber>::min();
		for (TInt ci = p_start ; ci < p_end; ++ci){
			if (vals[ci] >  maxinpartition[p]){
				maxinpartition[p] = vals[ci];
			}
		}
	}
}

/* same as inline void set_minexclusionnocheck in arrutilv2l1.h */
template <typename TInt, typename TNumber>
void set_min_excluding(TInt nvals, const TNumber * const vals, TInt excl, TNumber & toset){
	toset = std::numeric_limits<TNumber>::max();
	for (TInt i = 0; i < excl; ++i){
		if (vals[i] < toset){
			toset = vals[i];
		}
	}
	
	for (TInt i = excl+1; i < nvals; ++i){
		if (vals[i] < toset){
			toset = vals[i];
		}
	}
}

template <typename TSize, typename TInt>
TInt get_sum_int_array(TSize size, TInt * arr){
	TInt sum = arr[0];
	for (TSize i = 1; i < size; ++i){
		sum += arr[i];
	}
	return sum;
}

template <typename TSize, typename TInt>
bool get_sum_int_array_iszero(TSize size, const TInt * const arr){
	TInt sum = arr[0];
	TSize i = 0;
	while (sum == 0 && i < size){
		sum += arr[i];
		++i;
	}
	if (sum == 0){
		return true;
	}
	else{
		return false;
	}
}


template <typename TIntArray, typename TInt>
void integraladdto(TInt N, const TIntArray * const to_add, TIntArray * const to){
	for (TInt i = 0; i < N; ++i){
		to[i] += to_add[i];
	}
}


template <typename TInt>
std::vector<TInt> intlinspace(TInt i0, TInt i1, TInt npts){
	std::vector<TInt> linspaced (npts, i1);
	for (TInt k = 0; k < npts - 1; ++k){
		linspaced[k] = i0 + (k*(i1 - i0)/(npts - 1));
  }
  
  return linspaced;
}


template<typename TSize, typename TLabel>
void make_balanced(TSize minclustersize, TSize ndata, TLabel * const L, TSize nclusters, TSize * const groupsizes){
	
	TLabel argminc;
	TSize minc;
	TLabel argmaxc;
	TSize maxc;
	
	arrutilv2::set_argminmin(nclusters, groupsizes, argminc, minc);
	arrutilv2::set_argmaxmax(nclusters, groupsizes, argmaxc, maxc);
	
	while (minc < minclustersize){
		if (maxc <= minclustersize){
			throw std::runtime_error("In get contig by cluster 3, trying to balance clusters. minc < minclustersize, so balancing required. But maxc <= minclustersize, so balancing will cause another hole to appear, it will be impossible to fill all the holes there are simply not enough plugs");
		}
		
		*(std::find(L, L + ndata, argmaxc)) = argminc;
		++groupsizes[argminc];
		--groupsizes[argmaxc];
		arrutilv2::set_argminmin(nclusters, groupsizes, argminc, minc);
		arrutilv2::set_argmaxmax(nclusters, groupsizes, argmaxc, maxc);
	}
}


}

#endif
