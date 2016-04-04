/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef ARRUTILV2L1_H
#define ARRUTILV2L1_H

#include "arrutilv2l0.h"
#include "arrutilv2minmax.h"
#include <memory>

/*The functions below:
 * - would be the same if blas were used (not as low level as l0), sqrting, gramm matrices etc.
 * - should not use any kmeans words like centroids or data
 * */

namespace arrutilv2{	




template <typename TInt, typename TFloat>
inline void set_l2(const TInt & dimension, const TFloat * const a, const TFloat * const b, const TFloat & a_l22, const TFloat & b_l22, TFloat & l2){
	set_l22(dimension, a, b, a_l22, b_l22, l2);
	l2 = std::sqrt(std::max(static_cast<TFloat>(0), l2));
}



template <typename TInt, typename TFloat>
inline void set_l2(TInt dimension, const TFloat * const a, const TFloat * const b, const TFloat & a_l22, const TFloat & b_l22, TFloat & l2, TInt & ndcalcs){
	set_l2(dimension, a, b, a_l22, b_l22, l2);
	++ndcalcs;
}

template <typename TInt, typename TFloat>
inline TFloat get_l2(const TInt & dimension, const TFloat * const a, const TFloat * const b, const TFloat & a_l22, const TFloat & b_l22){
	TFloat l2;
	set_l2(dimension, a, b, a_l22, b_l22, l2);
	return l2;
}

template <typename TInt, typename TFloat>
inline TFloat get_l2(const TInt & dimension, const TFloat * const a, const TFloat * const b, const TFloat & a_l22, const TFloat & b_l22, TInt & ndcalcs){
	++ndcalcs;
	return get_l2(dimension, a, b, a_l22, b_l22);
}


/* input 
 * xxxxxxx (v, ncols) 
 * and
 * xxxxxxx (B, nrows x ncols)
 * xxxxxxx
 * xxxxxxx
 * output
 * xxx (nrows, distance from v to each row of B)
 */
template <typename TInt, typename TFloat>
void set_rl2s(TInt ncols, const TFloat * const v, TInt nrows, const TFloat * const B, const TFloat & v_l22, const TFloat * const B_l22s, TFloat * const l2s){
	set_rl22s(ncols,  v,  nrows,  B,  v_l22,  B_l22s, l2s);
	for(TInt r = 0; r < nrows; ++r){
		l2s[r] = std::sqrt(std::max(static_cast<TFloat>(0),l2s[r]));
	}
}




template <typename TInt, typename TFloat>
void set_rl2s(TInt ncols, const TFloat * const v, TInt nrows, const TFloat * const B, const TFloat & v_l22, const TFloat * const B_l22s, TFloat * const l2s, TInt & ndcalcs){
	set_rl2s(ncols, v, nrows,  B,  v_l22,  B_l22s, l2s);
	ndcalcs += nrows;
}

	
template <typename TInt, typename TFloat>
void set_rl22s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l22s){
	set_l22s(nrows, ncols, A, l22s, true);
}


template <typename TInt, typename TFloat>
void set_cl22s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l22s){
	set_l22s(nrows, ncols, A, l22s, false);
}




/* TODO : make byrow a template parameter */
template <typename TInt, typename TFloat>
inline void set_l2s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l2s, bool byrow){
	set_l22s(nrows, ncols, A, l2s, byrow);
	if (byrow == true){
		for (TInt r = 0; r < nrows; ++r){
			l2s[r] = std::sqrt(std::max(static_cast<TFloat>(0),l2s[r]));
		}
	}
	else{
		for(TInt c = 0; c < ncols; ++c){
			l2s[c] = std::sqrt(std::max(static_cast<TFloat>(0),l2s[c]));
		}
	}
}



template <typename TInt, typename TFloat>
void set_rl2s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l2s){	
	set_l2s(nrows, ncols, A, l2s, true);
}

template <typename TInt, typename TFloat>
void set_cl2s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l2s){	
	set_l2s(nrows, ncols, A, l2s, false);
}



template <typename TInt, typename TFloat>
void set_rsums(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const sums){	
	set_sums(nrows, ncols, A, sums, true);
}

template <typename TInt, typename TFloat>
void set_csums(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const sums){	
	set_sums(nrows, ncols, A, sums, false);
}

template <typename TInt, typename TFloat>
std::unique_ptr< TFloat []> getrsums(TInt nrows, TInt ncols, const TFloat * const A){
	std::unique_ptr< TFloat []> sums (new TFloat [nrows]);
	set_rsums(nrows, ncols, A, sums.get());
	return sums;
}


template <typename TInt, typename TFloat>
void set_l2s(TInt nrows, TInt ncols, TFloat * const old_array, const TFloat * const new_array, TFloat * const l2s, bool byrow){
	subtractfrom(nrows*ncols, new_array, old_array);
	if (byrow == true){
		set_rl2s(nrows, ncols, old_array, l2s);
	}
	
	else{
		set_cl2s(nrows, ncols, old_array, l2s);
	}
}



template <typename TInt, typename TFloat>
void set_rl2s(TInt nrows, TInt ncols, TFloat * const old_array, const TFloat * const new_array, TFloat * const l2s){
	set_l2s(nrows, ncols, old_array, new_array, l2s, true);
}

template <typename TInt, typename TFloat>
void set_cl2s(TInt nrows, TInt ncols, TFloat * const old_array, const TFloat * const new_array, TFloat * const l2s){
	set_l2s(nrows, ncols, old_array, new_array, l2s, false);
}



template <typename TInt, typename TFloat>
void set_l2sconst(TInt nrows, TInt ncols, const TFloat * const array1, const TFloat * const array2, TFloat * const l2s, bool byrow){

	//this memory allocation does not affect run time,
	std::unique_ptr< TFloat [] > array1_temp (new TFloat [nrows*ncols] );
	std::memcpy(array1_temp.get(), array1, nrows*ncols*sizeof(TFloat));
	set_l2s(nrows, ncols, array1_temp.get(), array2, l2s, byrow);
	
	
	
}




template <typename TInt, typename TFloat>
void set_rl2sconst(TInt nrows, TInt ncols, const TFloat * const array1, const TFloat * const array2, TFloat * const l2s){
	set_l2sconst(nrows, ncols, array1, array2, l2s, true);
}



template <typename TInt, typename TFloat>
void set_rl2sconst(TInt nrows, TInt ncols, const TFloat * const array1, const TFloat * const array2, TFloat * const l2s, TInt & ndcalcs){
	set_l2sconst(nrows, ncols, array1, array2, l2s, true);
	ndcalcs += nrows;
}


template <typename TInt, typename TFloat>
void set_cl2sconst(TInt nrows, TInt ncols, const TFloat * const array1, const TFloat * const array2, TFloat * const l2s){
	set_l2sconst(nrows, ncols, array1, array2, l2s, false);
}

/* C[i,j] <- distance squared from A[i,:] to B[j,:] */
template <typename TInt, typename TFloat>
void set_rrl22ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, TFloat * const l22ss){
	std::unique_ptr<TFloat []> A_ss (new TFloat [nrowsA]);	
	set_rl22s (nrowsA, ncols, A, A_ss.get());
	std::unique_ptr<TFloat []> B_ss (new TFloat [nrowsB]);	
	set_rl22s (nrowsB, ncols, B, B_ss.get()); //<TFloat, TInt>
	set_rrl22ss(nrowsA, ncols,  A,  nrowsB,  B, A_ss.get(), B_ss.get(),  l22ss);
}

template <typename TInt, typename TFloat>
void set_rrl2ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, TFloat * const l2ss){
	set_rrl22ss(nrowsA, ncols, A, nrowsB, B, l2ss);
	for (TInt i = 0; i < nrowsA*nrowsB; ++i){
		l2ss[i] = std::sqrt(std::max(static_cast<TFloat>(0), l2ss[i]));
	}
}

template <typename TInt, typename TFloat>
void set_rrl2ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l2ss){
	set_rrl22ss(nrowsA, ncols, A, nrowsB, B, A_l22s, B_l22s, l2ss);
	for (TInt i = 0; i < nrowsA*nrowsB; ++i){
		l2ss[i] = std::sqrt(std::max(static_cast<TFloat>(0), l2ss[i]));
	}
}

template <typename TInt, typename TFloat>
void set_rrl2ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l2ss, TInt & n_distance_calculations){	
	set_rrl2ss(nrowsA, ncols, A, nrowsB, B, A_l22s, B_l22s, l2ss);
	n_distance_calculations += nrowsA*nrowsB;
}

template <typename TInt, typename TFloat>
std::unique_ptr<TFloat []> getrrl2ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TInt & n_distance_calculations){
	std::unique_ptr<TFloat []> rrl2ss (new TFloat [nrowsA*nrowsB]);
	set_rrl2ss(nrowsA, ncols,  A, nrowsB,  B,  A_l22s,  B_l22s, rrl2ss.get(), n_distance_calculations);
	return rrl2ss;
}

template <typename TInt, typename TFloat>
void set_rrl2ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, TFloat * const l2ss, TInt & n_distance_calculations){
	set_rrl2ss(nrowsA, ncols, A, nrowsB, B, l2ss);
	n_distance_calculations += nrowsA*nrowsB;
}

/* for back compatibility (order of parameters breaks personal rules) */
template <typename TInt, typename TFloat>
void set_rrl22ss(TInt nrowsA, TInt nrowsB, TInt dimension, const TFloat * const B, const TFloat * const A, const TFloat * const B_l22s, const TFloat * const A_l22s, TFloat * const l22ss){
	set_rrl22ss(nrowsA, dimension, A, nrowsB, B, A_l22s, B_l22s, l22ss);
}

/* for back compatibility (order of parameters breaks personal rules) */
template <typename TInt, typename TFloat>
void set_rrl2ss(TInt nrowsA, TInt nrowsB, TInt dimension, const TFloat * const B, const TFloat * const A, const TFloat * const B_l22s, const TFloat * const A_l22s, TFloat * const l2ss){
	set_rrl22ss(nrowsA, nrowsB, dimension,  B,  A,  B_l22s,  A_l22s,  l2ss);	
	for (TInt i = 0; i < nrowsA; ++ i){
		for (TInt j = 0; j < nrowsB; ++ j){
			l2ss[i*nrowsB + j] = std::sqrt(std::max(static_cast<TFloat>(0), l2ss[i*nrowsB + j]));
		}
	}
}



/* untested */
template <typename TInt, typename TFloat>
void set_l2gramm(TInt nrows, TInt dimension, const TFloat * const C, const TFloat * const C_l22s, TFloat * const CC){
	set_rrl22ss(nrows, dimension, C, nrows, C, C_l22s, C_l22s, CC);
	
	for (TInt cr = 0; cr < nrows; ++cr){
		for (TInt cc = 0; cc < nrows; ++cc){
			CC[cr*nrows + cc] = std::sqrt(std::max(static_cast<TFloat>(0), CC[cr*nrows + cc]));

		}
	}
}

template <typename TInt, typename TFloat>
void set_l2gramm(TInt nrows, TInt dimension, const TFloat * const C, const TFloat * const C_l22s, TFloat * const CC, TInt & ndcalcs){
	set_l2gramm(nrows, dimension, C, C_l22s, CC);
	ndcalcs += nrows*nrows;
}



/* set_ rows r1-> r2 of the gramm matrix of CC (note : CC is pointer to full gramm matrix, variation of this may be needed, B?) */
template <typename TInt, typename TFloat>
void set_l2gramm_partial(TInt r1, TInt r2, TInt nrows, TInt dimension, const TFloat * const C, const TFloat * const C_l22s, TFloat * const CC, TInt & ndcalcs){
	set_rrl22ss(r2 - r1, dimension, C + r1*dimension, nrows, C, C_l22s + r1, C_l22s, CC + r1*nrows);
	for (TInt i = 0; i < (r2 - r1)*nrows; ++i){
		CC[r1*nrows + i] = std::sqrt(std::max(static_cast<TFloat>(0), CC[r1*nrows + i]));
	}
	ndcalcs += nrows*(r2 - r1);	
}



template <typename TInt, typename TFloat>
std::unique_ptr<TFloat []> get_rl22s(TInt nrows, TInt ncols, const TFloat * const A){
	std::unique_ptr<TFloat []> A_l22s (new TFloat[nrows]);
	set_rl22s(nrows, ncols, A, A_l22s.get());
	return A_l22s;	
}


template <typename TInt, typename TFloat>
std::unique_ptr<TFloat []> get_l2gramm(TInt nrows, TInt dimension, const TFloat * const C, const TFloat * const C_l22s, TInt & ndcalcs){
	std::unique_ptr<TFloat []> CC (new TFloat[nrows*nrows]); //Mr Meyers would disapprove for two reasons...
	set_l2gramm(nrows, dimension, C, C_l22s, CC.get(), ndcalcs);
	return CC;
}

template <typename TInt, typename TFloat>
std::unique_ptr<TFloat []> get_l2gramm(TInt nrows, TInt dimension, const TFloat * const C, const TFloat * const C_l22s){
	std::unique_ptr<TFloat []> CC (new TFloat[nrows*nrows]);
	set_l2gramm(nrows, dimension, C, C_l22s, CC.get());
	return CC;
}


/* argminss[i] = argmin_{ci \in nrowsB} | A[i] - B[ci] |  */
template <typename TInt, typename TFloat>
void set_rargmins(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TInt * const argmins){
	std::unique_ptr<TFloat []> l22s (new TFloat [nrowsB]);
	for (TInt i = 0; i < nrowsA; ++i){
		set_rl22s(ncols, A + i*ncols, nrowsB, B, A_l22s[i], B_l22s, l22s.get());
		set_argmin(nrowsB, l22s.get(),  argmins[i]);
	}	
}

template <typename TInt, typename TFloat>
void set_rargmins(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TInt * const argmins, TInt & ndcalcs){
	set_rargmins(nrowsA, ncols, A, nrowsB, B, A_l22s, B_l22s, argmins);
	ndcalcs += nrowsA*nrowsB;
}

template <typename TInt, typename TFloat>
void set_rargminmins(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TInt * const argmins, TFloat * const mins){
	std::unique_ptr<TFloat []> l22s (new TFloat [nrowsB]);
	for (TInt i = 0; i < nrowsA; ++i){
		set_rl22s(ncols, A + i*ncols, nrowsB, B, A_l22s[i], B_l22s, l22s.get());
		set_argminmin(nrowsB, l22s.get(),  argmins[i], mins[i]);
	}
}






template <typename TInt, typename TFloat>
void set_rargminmins(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TInt * const argmins, TFloat * const mins, TInt & ndcalcs){
	set_rargminmins(nrowsA, ncols, A, nrowsB, B, A_l22s, B_l22s, argmins, mins);
	ndcalcs += nrowsA*nrowsB;
}






}


#endif

//pushme


