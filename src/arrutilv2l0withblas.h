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

#ifndef ARRUTILV2L0WITHBLAS_H
#define ARRUTILV2L0WITHBLAS_H

#include "blastemplates.h"

/* blas versions of all functions in ../../blaslessbase/include/arrutilv2l0.h */

namespace arrutilv2{

inline void proxy_openblas_set_num_threads(int nthreads){
	openblas_set_num_threads(nthreads);
}



template <typename TInt, typename TFloat>
void set_rargmaxabs(TInt nrows, TInt ncols, const TFloat * const A,  TInt * argmaxabss){
	for (TInt r = 0; r < nrows; ++r){
		argmaxabss[r] = static_cast<TInt> (wblas::i_amax<TFloat>(ncols, A + r*ncols, 1));
	}
}


	
template <typename TInt, typename TFloat, typename TScaleFloatType>
void scale(TInt N, TScaleFloatType factor, TFloat * const toscale){
	TFloat factor_ct = static_cast<TFloat> (factor);
	wblas::scal(N, factor_ct, toscale, 1);
}



template <typename TInt, typename TFloat>
void rank1rowupdate(TInt ncols,  const TFloat * const row, TFloat scale, TInt nrows, TFloat * const toupdate){
	std::unique_ptr<TFloat []> ones (new TFloat [nrows]);//(nrows, 1);
	std::fill_n(ones.get(), nrows, 1); 
	wblas::ger(CblasRowMajor, nrows, ncols, scale, ones.get(), 1, row, 1, toupdate, ncols);
}		

template <typename TInt, typename TFloat>
inline void set_l22(const TInt & ndata, const TFloat * const a, TFloat & l22){
	/* for ndata = 28, slower, for ndata = 784, only slightly faster. Not going to use. */
	//l22 = wblas::dot<TFloat>(ndata, a, 1 , a, 1); 
	l22 = 0;
	for (TInt i = 0; i < ndata; ++i){
		l22 += a[i]*a[i];
	}
}

template <typename TInt, typename TFloat>
inline void set_l22(const TInt & dimension, const TFloat * const a, const TFloat * const b, const TFloat & a_l22, const TFloat & b_l22, TFloat & l22){
	/* see comment for previous set_l22 */
	l22 = 0;
	for (TInt i = 0; i < dimension; ++ i){
		l22 += a[i]*b[i];
	}
	l22 *= -2;
	l22 += a_l22;
	l22 += b_l22;
}	


template <typename TInt, typename TFloat>
inline void set_sum(const TInt & ndata, const TFloat * const a, TFloat & sum){
	/* unblasified, not worth it */
	sum = 0;
	for (TInt i = 0; i < ndata; ++i){
		sum += a[i];
	}
}


template <typename TInt, typename TFloat>
void set_l22s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l22s, bool byrow){
	/* tried blasifying, ran slower in d = 28 */
	if (byrow == true){
		for (TInt r = 0; r < nrows; ++r){
			set_l22(ncols, A + r*ncols, l22s[r]);
		}
	}
	
	else{
		for (TInt c = 0; c < ncols; ++c){
			l22s[c] = 0;
		}
		
		for (TInt r = 0; r < nrows; ++r){
			for (TInt c = 0; c < ncols; ++c){
				l22s[c] += A[r*ncols + c]*A[r*ncols + c];
			}
		}
	}
}




template <typename TInt, typename TFloat>
void set_sums(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const sums, bool byrow){
	/* TODO : blasify with a gemm */
	if (byrow == true){
		for (TInt r = 0; r < nrows; ++r){
			set_sum(ncols, A + r*ncols, sums[r]);
		}
	}
	
	else{
		for (TInt c = 0; c < ncols; ++c){
			sums[c] = 0;
		}
		for (TInt r = 0; r < nrows; ++r){
			for (TInt c = 0; c < ncols; ++c){
				sums[c] += A[r*ncols + c];
			}
		}
	}
}

/* v : 1-D of size ncols
* B : nrows x ncols
* l22s[i] : |v - B[i]|_2
* */
template <typename TInt, typename TFloat>
inline void set_rl22s(const TInt & ncols, const TFloat * const v, const TInt & nrows, const TFloat * const B, const TFloat & v_l22s, const TFloat * const B_l22s, TFloat * const l22s){	
	//-2 cross term. using blas really helps here
	wblas::gemv<TFloat>(CblasRowMajor, CblasNoTrans, nrows, ncols, -2, B, ncols, v, 1, 0, l22s, 1);
	for (TInt r = 0; r < nrows; ++r){
		l22s[r] += v_l22s;
		l22s[r] += B_l22s[r];
	}
}


/* A : nrowsA x ncols
* B : nrowsB x ncols
* C[i,j] : |A[i] - B[j]|_2 for 0 <= i < nrowsA   0 <= j < nrowsB 
* untested 
* */
template <typename TInt, typename TFloat>
void set_rrl22ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l22ss){
	wblas::gemm<TFloat>(CblasRowMajor, CblasNoTrans, CblasTrans, nrowsA, nrowsB, ncols, -2.0, A, ncols, B, ncols, 0.0, l22ss, nrowsB);
	TInt mlen = (nrowsA > nrowsB) ? nrowsA : nrowsB;
	std::unique_ptr<TFloat []> ones (new TFloat[mlen]);
	for (TInt i = 0; i < mlen; ++i){
		ones[i] = 1;
	}
	wblas::ger<TFloat>(CblasRowMajor, nrowsA, nrowsB, 1., ones.get(), 1, B_l22s, 1, l22ss, nrowsB);
	wblas::ger<TFloat>(CblasRowMajor, nrowsA, nrowsB, 1., A_l22s, 1, ones.get(), 1, l22ss, nrowsB);
}


template <typename TInt, typename TFloat>
void subtractfrom(TInt N, const TFloat * const tosubtract, TFloat * const from){
	/* with d = 28, blas version is slower. with d = 784, no difference. not going to use blas */
	//wblas::axpy<TFloat>(N, -1, tosubtract, 1, from, 1);
	for (TInt i = 0; i < N; ++i){
		from[i] -= tosubtract[i];
	}
}	

template <typename TInt, typename TFloat>
void addto(TInt N, const TFloat * const toadd, TFloat * const to){
	/* see subtractfrom comment */
	for (TInt i = 0; i < N; ++i){
		to[i] += toadd[i];
	}
}






}

#endif




