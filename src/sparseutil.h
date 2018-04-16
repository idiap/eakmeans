/*
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <james.newling@gmail.com>
All rights reserved.

eakmeans is a library for exact and approximate k-means written in C++ and
Python. This file is part of eakmeans. See file COPYING for more details.

This file is part of eakmeans.

eakmeans is free software: you can redistribute it and/or modify
it under the terms of the 3-Clause BSD Licence. See
https://opensource.org/licenses/BSD-3-Clause for more details.

eakmeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
COPYING for more details.
*/

#ifndef SPARSEUTIL_H
#define SPARSEUTIL_H

#include <algorithm>

#include "sparsedatasets.h"
#include <vector>

#include <cmath>

#include <memory>
#include <mutex>
#include <atomic>

#include <algorithm>
#include <iomanip>
/* 
 * */
 
namespace sparse{

template <typename TInt, typename TFloat>
inline void set_label(TInt n_nonzero, const TInt * const indices, const TFloat * const values, TInt dimension, TInt ncentroids, const TFloat * const C, const TFloat * const C_l22s, TInt & label){ //, const TFloat & data_l22 not nec, just a constant factor.

	TFloat nearest_uptoc = 10e44;
	TFloat distance2_uptoc;
	
	for (TInt ci = 0; ci < ncentroids; ++ ci){
		distance2_uptoc = 0;		
		for (TInt di = 0; di < n_nonzero; ++di){
			distance2_uptoc += values[di]*C[ci*dimension + indices[di]];
		}
		distance2_uptoc *= -2.;
		//distance2 += data_l22;//not nec, same for all.
		distance2_uptoc += C_l22s[ci];
				
		if (distance2_uptoc < nearest_uptoc){
			label = ci;
			nearest_uptoc = distance2_uptoc;
		}
	}
}




/* Given sparse vector (n_nonzero, indices, values) and centroids, set argmin, min and l2s. */
template <typename TInt, typename TFloat>
inline void set_argminmin_rl2s(TInt n_nonzero, const TInt * const indices, const TFloat * const values, 
TInt dimension, TInt ncentroids, const TFloat * const C, 
const TFloat & data_l22, const TFloat * const C_l22s, 
TInt & argmin, TFloat & min, TFloat * const l2s){

	min = 10e44;
	TFloat distance2;
	
	for (TInt ci = 0; ci < ncentroids; ++ ci){
		distance2 = 0;		
		for (TInt di = 0; di < n_nonzero; ++di){
			distance2 += values[di]*C[ci*dimension + indices[di]];
		}
		distance2 *= -2.;
		distance2 += data_l22;
		distance2 += C_l22s[ci];
		
		l2s[ci] = std::sqrt(std::max(static_cast<TFloat>(0), distance2));
		if (l2s[ci] < min){
			argmin = ci;
			min = l2s[ci];
		}
	}
}


 
//set L in range [data0, data1).
template <typename TInt, typename TFloat>
void update_L(const SparseData<TInt, TFloat> & sd, 
TInt data0, TInt data1, 
TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, 
TInt * const L, std::vector<std::tuple<TInt, TInt, TInt> > & where_label_change){
	TInt old_label;
	
	const TInt * const starts =  sd.starts.data();
	const TInt * const indices =  sd.indices.data();
	const TFloat * const values =  sd.values.data();
	const TInt dimension = sd.dimension;
	
	TFloat nearest_uptoc;
	TFloat distance2_uptoc;
	
	for (TInt i = data0; i < data1; ++i){
		old_label = L[i];
	
		/* what follows is essentially inlining of set_label, this forced inlining accelerates running time by 30%.
		 * the original set_label call was 
		 * set_label(starts[i+1] - starts[i], indices + starts[i], values + starts[i], dimension, ncentroids, C,  C_l22s, L[i]);
		 * */
		nearest_uptoc = 10e44;
		for (TInt ci = 0; ci < ncentroids; ++ ci){
			distance2_uptoc = 0;		
			for (TInt di = starts[i]; di < starts[i+1]; ++di){
				distance2_uptoc += values[di]*C[ci*dimension + indices[di]];
			}
			distance2_uptoc *= -2.;
			distance2_uptoc += C_l22s[ci];
				
			if (distance2_uptoc < nearest_uptoc){
				L[i] = ci;
				nearest_uptoc = distance2_uptoc;
			}
		}
				
		if (L[i] != old_label){
			where_label_change.emplace_back(i, old_label, L[i]);
		}
	}
}

//update L and d in range [data0, data1), no unput variables are offset.
template <typename TInt, typename TFloat>
void update_L_dn(const SparseData<TInt, TFloat> & sd, 
TInt data0, TInt data1, 
TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, 
TInt * const L, TFloat * const dn, std::vector<std::tuple<TInt, TInt, TInt> > & where_label_change){
		
	TInt old_label;
	const TInt * const starts =  sd.starts.data();
	const TInt * const indices =  sd.indices.data();
	const TFloat * const values =  sd.values.data();
	const TInt dimension = sd.dimension;

	//distance squared from a centroid to the datapoint, minus the norm of the datapoint squared (constant)
	TFloat distance2_uptoconstant;
	//the lowest over centroids of above value found so far. 
	TFloat d2nearest_uptoconstant;

	
	for (TInt i = data0; i < data1; ++i){
		old_label = L[i];
		d2nearest_uptoconstant = 10e44;
		for (TInt ci = 0; ci < ncentroids; ++ ci){
			distance2_uptoconstant = 0;		
			for (TInt di = starts[i]; di < starts[i+1]; ++di){
				distance2_uptoconstant += values[di]*C[ci*dimension + indices[di]];
			}
			distance2_uptoconstant *= -2.;
			distance2_uptoconstant += C_l22s[ci];
				
			if (distance2_uptoconstant < d2nearest_uptoconstant){
				L[i] = ci;
				d2nearest_uptoconstant = distance2_uptoconstant;
			}
		}
		
		//add the missing constant and square root it
		dn[i] = d2nearest_uptoconstant + data_l22s[i];
		dn[i] = std::sqrt(std::max(static_cast<TFloat>(0), dn[i])); 
		
		if (L[i] != old_label){
			where_label_change.emplace_back(i, old_label, L[i]);
		}	
	}
}


//where no variables are offset, set L and dn in the range (data0, data1)
template <typename TInt, typename TFloat>
void set_L_dn(const SparseData<TInt, TFloat> & sd, 
TInt data0, TInt data1, 
TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, 
TInt * const L, TFloat * const dn){


	const TInt * const starts =  sd.starts.data();
	const TInt * const indices =  sd.indices.data();
	const TFloat * const values =  sd.values.data();
	const TInt dimension = sd.dimension;
	
	//distance squared from a centroid to the datapoint, minus the norm of the datapoint squared (constant)
	TFloat distance2_uptoconstant;
	//the lowest over centroids of above value found so far. 
	TFloat d2nearest_uptoconstant;

	
	for (TInt i = data0; i < data1; ++i){
		
		d2nearest_uptoconstant = 10e44;
		for (TInt ci = 0; ci < ncentroids; ++ ci){
			distance2_uptoconstant = 0;		
			for (TInt di = starts[i]; di < starts[i+1]; ++di){
				distance2_uptoconstant += values[di]*C[ci*dimension + indices[di]];
			}
			distance2_uptoconstant *= -2.;
			distance2_uptoconstant += C_l22s[ci];
				
			if (distance2_uptoconstant < d2nearest_uptoconstant){
				L[i] = ci;
				d2nearest_uptoconstant = distance2_uptoconstant;
			}
		}
		
		//add the missing constant,
		dn[i] = d2nearest_uptoconstant + data_l22s[i]; 
	}
}



//where no variables are offset, set L, lowers (all distances) and dn in the range (data0, data1)
template <typename TInt, typename TFloat>
void set_L_lowers_dn(const SparseData<TInt, TFloat> & sd, 
TInt data0, TInt data1, 
TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, 
TInt * const L, TFloat * const lowers, TFloat * const dn){


	const TInt * const starts =  sd.starts.data();
	const TInt * const indices =  sd.indices.data();
	const TFloat * const values =  sd.values.data();
	const TInt dimension = sd.dimension;
	
	//distance squared from a centroid to the
	TFloat distance2;
	//the lowest over centroids of above value found so far. 
	TFloat d2nearest;

	
	for (TInt i = data0; i < data1; ++i){
		
		d2nearest = 10e44;
		for (TInt ci = 0; ci < ncentroids; ++ ci){
			distance2 = 0;		
			for (TInt di = starts[i]; di < starts[i+1]; ++di){
				distance2 += values[di]*C[ci*dimension + indices[di]];
			}
			distance2 *= -2.;
			distance2 += C_l22s[ci];
			distance2 += data_l22s[i];
			
			lowers[i*ncentroids + ci] = std::sqrt(std::max(static_cast<TFloat>(0), distance2));	
			if (distance2 < d2nearest){
				L[i] = ci;
				d2nearest = distance2;
			}			
		}
	}
}






template <typename TInt, typename TFloat>
void update_S_H(const SparseData<TInt, TFloat> & sd, const std::vector<std::tuple<TInt, TInt, TInt> > & where_label_change, TFloat * const S, TInt * const H){
	TInt where;
	for (auto & x: where_label_change){
		where =std::get<0>(x);
		for (TInt j = sd.starts[where]; j < sd.starts[where + 1]; ++j){
			S[sd.dimension*std::get<1>(x) + sd.indices[j]] -= sd.values[j];
			S[sd.dimension*std::get<2>(x) + sd.indices[j]] += sd.values[j];
		}

		--H[std::get<1>(x)];
		++H[std::get<2>(x)];
	}
	
}






template <typename TInt, typename TFloat>
void increment_S_H(TInt data0, TInt data1, const SparseData<TInt, TFloat> & sd, const TInt * const L, TFloat * const S, TInt * const H){
	
	for (TInt i = data0; i < data1; ++i){
		for (TInt j = sd.starts[i]; j < sd.starts[i + 1]; ++j){
			S[sd.dimension*L[i] + sd.indices[j]] += sd.values[j];
		}
		++H[L[i]];
	}
}



template <typename TInt, typename TFloat>
std::function<void(TInt)> update_L_label_changes_ati(TInt nthreads, const SparseData<TInt, TFloat> & data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, std::atomic<TInt> & ndcalcs, std::vector<std::vector<std::tuple<TInt, TInt, TInt>>> & where_label_changes){
	
	
	
		return [nthreads, &data, ncentroids, C, data_l22s, C_l22s, L, &ndcalcs, &where_label_changes](TInt ti){
			
			TInt x0 = (ti*data.ndata)/nthreads;
			TInt x1 = ((ti+1)*data.ndata)/nthreads;
			
			where_label_changes[ti].clear(); //index, old, new.
			update_L(data, x0, x1, ncentroids, C, data_l22s, C_l22s, L, where_label_changes[ti]);
			ndcalcs += ncentroids*(x1 - x0);
		
		};
	//}	
	
}


template <typename TInt, typename TFloat>
void update_S_H_from_label_changes(
const SparseData<TInt, TFloat> & data, std::vector<std::vector<std::tuple<TInt, TInt, TInt>>> & where_label_changes, TFloat * const sums, TInt * const counts){
	for (auto & x: where_label_changes){
		update_S_H(data, x, sums, counts);
	}
}


template <typename TInt, typename TFloat>
void update_S_H_from_label_changes(
const SparseData<TInt, TFloat> & data, std::vector<std::vector<std::tuple<TInt, TInt, TInt>>> & where_label_changes, TFloat * const sums, TInt * const counts, TInt & nchanges){
	for (auto & x: where_label_changes){
		update_S_H(data, x, sums, counts);
		nchanges += x.size();
	}
}


	

template <typename TInt, typename TFloat>
std::function<void(TInt)> update_L_S_H_ati(TInt nthreads, const SparseData<TInt, TFloat> & sd, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const S, TInt * const H, TInt & nchanges, std::mutex & work_mutex, std::atomic<TInt> & ndcalcs, std::vector<std::vector<std::tuple<TInt, TInt, TInt>>> & where_label_changes)
{
	
	if (nthreads != 1){
		throw std::runtime_error("Currently only implemented for 1 thread " );
	}
	
	else{
		return [nthreads, &sd, ncentroids, C, data_l22s, C_l22s, L, S, H, &nchanges, &work_mutex, &ndcalcs, &where_label_changes](TInt){
			
			//
			where_label_changes[0].clear(); //index, old, new.
			update_L(sd, static_cast<TInt> (0), sd.ndata, ncentroids, C, data_l22s, C_l22s, L, where_label_changes[0]);
			ndcalcs += ncentroids*sd.ndata;
			nchanges += where_label_changes[0].size();
			
			update_S_H(sd, where_label_changes[0], S, H);

		};
	}	
}



template <typename TInt, typename TFloat> 
void set_inner(
TInt n_a, const  TInt * const a_indices, const TFloat *  const a_values,  
TInt n_b, const TInt * const b_indices, const TFloat *  const b_values, 
TFloat & inner){
	TInt a_i = 0;
	TInt b_i = 0;
	inner = 0;
	while (a_i < n_a && b_i < n_b){
		
		if (a_indices[a_i] < b_indices[b_i]){
			++a_i;
		}
		else if (b_indices[b_i] < a_indices[a_i]){
			++b_i;
		}
		else{
			inner += a_values[a_i]*b_values[b_i];
			++a_i;
			++b_i;
		}
		
	}
}

template <typename TInt, typename TFloat> 
TFloat get_inner(
TInt n_a, const TInt * const  a_indices, const TFloat *  const a_values,  
TInt n_b, const TInt * const  b_indices, const TFloat *  const b_values){
	TFloat inner;
	set_inner(n_a, a_indices, a_values, n_b, b_indices, b_values, inner);
	return inner;	
}


	

template <typename TInt, typename TFloat> 
void set_l22(
TInt n_a, const TFloat *  const a_values, TFloat & l22){ //const TInt * const  a_indices, 
	l22 = 0;
	for (TInt j = 0; j < n_a; ++j){
		l22 += a_values[j]*a_values[j];
	}
}

template <typename TInt, typename TFloat> 
TFloat get_l22(
TInt n_a, const TInt * const  a_indices, const TFloat *  const a_values,
TInt n_b, const TInt * const  b_indices, const TFloat *  const b_values)
{
	TFloat a_l22, b_l22;
	set_l22(n_a,  a_values, a_l22); //a_indices,  
	set_l22(n_b,  b_values, b_l22); // b_indices,  
	return a_l22 + b_l22 - 2*get_inner(n_a, a_indices, a_values, n_b, b_indices, b_values);
}




template <typename TInt, typename TFloat> 
void set_rl22s(const sparse::SparseData<TInt, TFloat>  & sd, TFloat * const rl22s){
	for (TInt i = 0; i < sd.ndata; ++i){
		rl22s[i] = 0;
		for (TInt j = sd.starts[i]; j < sd.starts[i+1]; ++j){
			rl22s[i] += sd.values[j]*sd.values[j];
		}
	}
}

//inner between sparse and dense
template <typename TInt, typename TFloat>
inline TFloat get_inner( 
TInt n_a, const TInt * const  a_indices, const TFloat *  const a_values,
const TFloat * const b 
){
	TFloat inner = 0;
	for (TInt j = 0; j < n_a; ++j){
		inner += a_values[j]*b[a_indices[j]];
	}
	return inner;
}



template <typename TInt, typename TFloat> 
std::unique_ptr<TFloat []> get_rl22s(const sparse::SparseData<TInt, TFloat>  & sd){
	std::unique_ptr<TFloat []> rl22s (new TFloat [sd.ndata]);
	set_rl22s(sd, rl22s.get());
	return rl22s;
}






template <typename TInt, typename TFloat>
TFloat getmeanl22at(const SparseData<TInt, TFloat> & sd, const TFloat * const C, const TInt * const L, const TFloat * const data_l22s, const TFloat * const C_l22s){
	TFloat sse = 0;
	for (TInt i = 0; i < sd.ndata; ++i){
		
		
		sse += 
		data_l22s[i] + C_l22s[L[i]]
		-2.*get_inner(sd.starts[i+1] - sd.starts[i], 
		sd.indices.data() + sd.starts[i], 
		sd.values.data() + sd.starts[i],
		C + sd.dimension*L[i]);
	}
	TFloat mse = sse / static_cast<TFloat> (sd.ndata); 
	
	return  mse;
}
	



template <typename TInt, typename TFloat>
inline void set_l22(TInt n_nonzero, const TInt * const a_indices, const TFloat * const a_values, const TFloat * const b_dense, const TFloat & a_l22, const TFloat & b_l22, TFloat & l22, TInt & ndcalcs){
	ndcalcs += 1;
	l22 = a_l22 + b_l22 - 2*get_inner(n_nonzero, a_indices, a_values, b_dense);	
}


//template <typename TInt, typename TFloat>
//inline void set_l2(TInt n_nonzero, const TInt * const a_indices, const TFloat * const a_values, const TFloat * const b_dense, const TFloat & a_l22, const TFloat & b_l22, TFloat & l2, TInt & ndcalcs){
	//ndcalcs += 1;
	//l2 = std::sqrt(
	//std::max<TFloat>(static_cast<TFloat>(0), 
	//a_l22 + b_l22 - 2*get_inner(n_nonzero, a_indices, a_values, b_dense))
	//);
//}

template <typename TInt, typename TFloat>
inline void set_l2(TInt n_nonzero, const TInt * const a_indices, const TFloat * const a_values, const TFloat * const b_dense, const TFloat & a_l22, const TFloat & b_l22, TFloat & l2, TInt & ndcalcs){
	
	TFloat inner = 0;//a_values[0]*b_dense[a_indices[0]];
	for (TInt j = 0; j < n_nonzero; ++j){
		inner += a_values[j]*b_dense[a_indices[j]];
	}	
	l2 = std::sqrt(std::max<TFloat>(static_cast<TFloat>(0), a_l22 + b_l22 - 2*inner));
	ndcalcs += 1;

}


	  



namespace todense{

template <typename TInt, typename TFloat>
void zero_and_copy(TInt index, const sparse::SparseData<TInt, TFloat> & sd, TFloat * const a){
	std::fill_n(a, sd.dimension, 0); 
	for (TInt j = sd.starts[index]; j < sd.starts[index +1]; ++j){
		a[sd.indices[j]] = sd.values[j];
	}
	
}
 
template <typename TInt, typename TFloat>
void increment_S_H(const SparseData<TInt, TFloat> & sd, TInt data0, TInt data1, const TInt * const L,  TFloat * const S,  TInt * const H){
	TInt dimension = sd.dimension;
	for (TInt i = data0; i < data1; ++i){
		TInt label = L[i];
		++H[label];
		for (TInt j = sd.starts[i]; j < sd.starts[i+1]; ++j){
			S[label*dimension + sd.indices[j]] += sd.values[j];
		}
	}	
}


template <typename TInt, typename TFloat>
void set_S_H(const SparseData<TInt, TFloat> & sd, TInt data0, TInt data1, TInt ncentroids, 
const TInt * const L,  TFloat * const S,  TInt * const H){


	TInt dimension = sd.dimension;
		
	std::fill_n(H, ncentroids, 0);
	std::fill_n(S, ncentroids*dimension, 0);

	
	TInt label;

	for (TInt i = data0; i < data1; ++i){
		label = L[i];
		++H[label];
		
	
		for (TInt j = sd.starts[i]; j < sd.starts[i+1]; ++j){
			S[label*dimension + sd.indices[j]] += sd.values[j];
		}
	}
	
	
}

template <typename TInt, typename TFloat>
void copyatindices(TInt n_A, const SparseData<TInt, TFloat> & sd, TFloat * const A, const TInt * const indices){
	
	TInt dimension = sd.dimension;
	std::fill_n(A, dimension*n_A, 0); // is this nec?
	for (TInt ci = 0; ci < n_A; ++ci){
		for (TInt j = sd.starts[indices[ci]]; j <  sd.starts[indices[ci] +1]; ++j){
			A[ci*dimension + sd.indices[j]] = sd.values[j];
		}		
	}
}



}
}


#endif
