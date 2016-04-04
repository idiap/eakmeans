/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef ARRUTILV2COPY_H
#define ARRUTILV2COPY_H

#include <exception>
#include <stdexcept>
#include <cstring>
#include <vector>
#include <numeric>
#include <memory>


namespace arrutilv2{
	
template <typename TSize, typename DType>
void copyatindices(TSize ndata_from, TSize ndata_to, TSize dimension, const DType * const data_from, DType * const data_to, const TSize * const indices){
		
	for (TSize index_i = 0; index_i < ndata_to; ++index_i){
		std::memcpy(data_to + index_i*dimension, data_from + indices[index_i]*dimension, dimension*sizeof(DType));
	}
}


/* calls copy at indices but first does some checks for uniqueness, bounds etc*/
template <typename TSize, typename DType>
void copyatuniqueindices(TSize ndata_from, TSize ndata_to, TSize dimension, const DType * const data_from, DType * const data_to, const TSize * const indices){
	
	std::vector<bool> used (ndata_from,false);
	TSize index;
	for (TSize i = 0; i < ndata_to; ++i){
		index = indices[i];
		if (index >= ndata_from){
			throw std::runtime_error("index of the data to copy (" + std::to_string(index) + ") seems to be out of range in copyatuniqueindices (" + std::to_string(ndata_from) + "). In other words, the following is false, causing the throw " + std::to_string(index) + " >= " + std::to_string(ndata_from));
		}
		
		else if (used[index] == true){
			throw std::runtime_error("index to copy (" +std::to_string(index) +") seems to have already been copied. Possible cause: indices to be copied are not unique, which contradicts assumptions in this function");
		}
		
		else{
			used[i] = true;
		}
	}
	copyatindices(ndata_from, ndata_to, dimension, data_from, data_to, indices);
}
	


template <typename TSize, typename DType>
void copyatindices(TSize ndata_from, TSize dimension, const DType * const data_from, DType * const data_to, const std::vector<TSize> & indices){
	TSize ndata_to = indices.size();
	copyatindices(ndata_from, ndata_to, dimension, data_from, data_to, indices.data());
}


template <typename TSize, typename DType>
void copyatuniqueindices(TSize ndata_from, TSize dimension, const DType * const data_from, DType * const data_to, const std::vector<TSize> & indices){
	TSize ndata_to = indices.size();
	copyatuniqueindices(ndata_from, ndata_to, dimension, data_from, data_to, indices.data());
}

//untested
template <typename TSize, typename DType>
std::unique_ptr<DType [] > getatuniqueindices(TSize ndata_from, TSize dimension, const DType * const data_from, const std::vector<TSize> & indices){
	std::unique_ptr<DType [] > data_to (new DType[dimension*indices.size()]);
	copyatuniqueindices(ndata_from, dimension, data_from, data_to.get(), indices);
	return data_to;
}

template <typename TSize, class T>
std::unique_ptr<T []> copy_uptrarr_to_uptrarr(TSize n, const std::unique_ptr<T []> & uptr){
	std::unique_ptr<T []> thecopy (new T[n]);
	std::memcpy(thecopy.get(), uptr.get(), n*sizeof(T));
	return thecopy;
}


template <typename TSize, class T>
std::unique_ptr<T []> copy_ptrarr_to_uptrarr(TSize n, const T * const ptrarr){
	std::unique_ptr<T []> thecopy (new T[n]);
	std::memcpy(thecopy.get(), ptrarr, n*sizeof(T));
	return thecopy;
}


template <typename TNumber, typename TSize>
std::unique_ptr<TNumber []> get_initialised_uptrarr(TSize npts, TNumber defval){
	std::unique_ptr<TNumber []> upined (new TNumber [npts]);
	std::fill_n(upined.get(), npts, defval);
	return upined; 
}


template <typename TInt>
std::unique_ptr < TInt []> get_with_offset(TInt n, const TInt * const vals, TInt offset){
	std::unique_ptr < TInt []> newvals (new TInt [n]);
	for (TInt ci = 0; ci < n; ++ ci){
		newvals[ci] = offset + vals[ci];
	}
	return newvals;
}


}

#endif
