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

#ifndef ARRUTILV2_MINMAX_H
#define ARRUTILV2_MINMAX_H

#include <memory>
#include <numeric> 

/*The functions below:
 * - low level mins and max, independent of everything else. 
 * - set one variable or argmin&min*min2 (no multiple settings)
 *  */

namespace arrutilv2{


template <typename TInt, typename TFloat>
TFloat getmaxval(TInt N, const TFloat * const a){
	TFloat mv = a[0];
	for (TInt i = 1; i < N; ++i){
		if (a[i] > mv){
			mv = a[i];
		}
	}
	return mv;
}

template <typename TInt, typename TFloat>
inline void set_min(const TInt & dimension, const TFloat * const a, TFloat & min){
	min = a[0];
	for (TInt i = 1; i < dimension; ++i){
		if (a[i] < min){
			min = a[i];
		}
	}
}


template <typename TInt, typename TFloat>
inline void set_minexclusionnocheck(const TInt & dimension, const TFloat * const a, const TInt & exclude, TFloat & min){
			
	min = std::numeric_limits<TFloat>::max();
	for (TInt i = 0; i < exclude; ++i){
		if (a[i] < min){
			min = a[i];
		}
	}
	for (TInt i = exclude + 1; i < dimension; ++i){
		if (a[i] < min){
			min = a[i];
		} 
	}
}


template <typename TInt, typename TFloat>
inline void set_minexclusion(const TInt & dimension, const TFloat * const a, const TInt & exclude, TFloat & min){
	if (dimension < 2){
		throw std::logic_error("cannot set min with exclusion if vector has length less than 2 in set_minexclusion");
	}
	
	if (exclude >= dimension){
		throw std::logic_error("exclusion index exceeds range in set_minexlusion");
	}
	
	set_minexclusionnocheck(dimension, a, exclude, min);
}	



template <typename TInt, typename TFloat>
inline void set_argminmin(const TInt & dimension, const TFloat * const a,  TInt & argmin, TFloat & min){
	argmin = 0;
	min = a[0];
	for (TInt i = 1; i < dimension; ++i){
		if (a[i] < min){
			min = a[i];
			argmin = i;
		}
	}
}

template <typename TInt, typename TFloat>
inline void set_argmaxmax(const TInt & dimension, const TFloat * const a,  TInt & argmax, TFloat & max){
	argmax = 0;
	max = a[0];
	for (TInt i = 1; i < dimension; ++i){
		if (a[i] > max){
			max = a[i];
			argmax = i;
		}
	}
}


template <typename TInt, typename TFloat>
inline void set_argmin(const TInt & dimension, const TFloat * const a,  TInt & argmin){
	TFloat min;
	set_argminmin(dimension, a, argmin, min);
}

template <typename TInt, typename TFloat>
inline TInt get_argmin(const TInt & dimension, const TFloat * const a){
	TFloat min;
	TInt argmin;
	set_argminmin(dimension, a, argmin, min);
	return argmin;
}




template <typename TInt, typename TFloat>
void set_rargminmins(TInt nrows, TInt ncols, const TFloat * const A,  TInt * const argmins, TFloat * const mins){
	for (TInt r = 0; r < nrows; ++r){
		set_argminmin(ncols, A + r*ncols, argmins[r], mins[r]);
	}
}

template <typename TInt, typename TFloat>
void set_rargmins(TInt nrows, TInt ncols, const TFloat * const A,  TInt * const argmins){
	for (TInt r = 0; r < nrows; ++r){
		set_argmin(ncols, A + r*ncols, argmins[r]);
	}
}



template <typename TInt, typename TFloat>
inline void set_argminmin2nocheck(const TInt & n, const TFloat * const a,  TInt & argmin, TFloat & min, TFloat & min2){

	
	if (a[0] < a[1]){
		min = a[0];
		argmin = 0;
		min2 = a[1];
	}
	
	else{
		min = a[1];
		argmin = 1;
		min2 = a[0];
	}
	
	for (TInt i = 2; i < n; ++i){
		if (a[i] < min){
			min2  = min;
			argmin = i;
			min = a[i];
		}
		
		else if (a[i] < min2){
			min2 = a[i];
		}
	}
}


template <typename TInt, typename TFloat>
inline void set_argmin2min2nocheck(const TInt & n, const TFloat * const a,  TInt & argmin, TInt & argmin2, TFloat & min, TFloat & min2){

	
	if (a[0] < a[1]){
		min = a[0];
		argmin = 0;
		min2 = a[1];
		argmin2 = 1;
	}
	
	else{
		min = a[1];
		argmin = 1;
		min2 = a[0];
		argmin2 = 0;
	}
	
	for (TInt i = 2; i < n; ++i){
		if (a[i] < min){
			
			argmin2 = argmin;
			min2  = min;
			argmin = i;
			min = a[i];
		}
		
		else if (a[i] < min2){
			argmin2 = i;
			min2 = a[i];
		}
	}
}


template <typename TInt, typename TFloat>
inline void set_argminmin2(const TInt & n, const TFloat * const a,  TInt & argmin, TFloat & min, TFloat & min2){
	if (n == 0 || n == 1){
		throw std::logic_error("attempt to set min and argmin and min2 in set_argminmin2, but the vector is of length 0 or 1.");
	}
	set_argminmin2nocheck(n, a, argmin, min, min2);
}


template <typename TInt, typename TFloat>
void set_rminsexclusion(TInt nrows, TInt ncols, const TFloat * const A, const TInt * const exclusions, TFloat * const min){
	if (ncols < 2){
		throw std::logic_error("ncols less than 2 in rminsexclusion");
	}
	
	if (*(std::max_element(exclusions, exclusions + nrows)) >= ncols){
		throw std::logic_error("max exclusion element in rminsexclusion is greater than or equal to ncols");
	}
		
	for (TInt i = 0; i < nrows; ++i){
		set_minexclusionnocheck(ncols, A + i*ncols, exclusions[i], min[i]);
	}
}	

template <typename TInt, typename TFloat>
void set_rargminminsexclusion(TInt nrows, TInt ncols, const TFloat * const A, const TInt * const exclusions, TInt * const argmin, TFloat * const min){
	for (TInt r = 0; r < nrows; ++r){
		if (exclusions[r] >= ncols){
			throw std::runtime_error("exclusion index exceeds number of columns in set__row_argminmin_exclusion");
		}
		min[r] = *(A + r*ncols) + *(A + r*ncols + 1) + 1;
		argmin[r] = 999;		
		for (TInt c = 0; c < exclusions[r]; ++ c){
			if (*(A + r*ncols + c) < min[r]){
				min[r] = *(A + r*ncols + c);
				argmin[r] = c;
			}
		}
		for (TInt c = exclusions[r] + 1; c < ncols; ++ c){
			if (*(A + r*ncols + c) < min[r]){
				min[r] = *(A + r*ncols + c);
				argmin[r] = c;
			}
		}
	}
}







template <typename TInt, typename TFloat>
void set_rargminmin2s(TInt nrows, TInt ncols, const TFloat * const A,  TInt * const argmin, TFloat * const min, TFloat * const secondmin){
	if (ncols == 0 || ncols == 1){
		throw std::runtime_error("call to set_rargminmin2s with ncols = 1 or 0 not allowed (at least 2 cols needed)");
	}
		
	else{
		for (TInt r = 0; r < nrows; ++r){
			set_argminmin2nocheck(ncols, A + r*ncols, argmin[r], min[r], secondmin[r]);
		}
	}
}





template <typename TInt, typename TFloat>
void set_rargminminsnodiag(TInt npts, const TFloat * const A,  TInt * const argmin, TFloat * const min){
	std::unique_ptr<TInt []> excl(new TInt[npts]);
	std::iota(excl.get(), excl.get() + npts, 0);
	set_rargminminsexclusion(npts, npts, A, excl.get(), argmin, min);
}
	



template <typename TInt, typename TFloat>
void set_rminsnodiag(TInt npts, const TFloat * const A, TFloat * const min){
	std::unique_ptr<TInt []> excl(new TInt[npts]);
	std::iota(excl.get(), excl.get() + npts, 0);
	set_rminsexclusion(npts, npts, A, excl.get(), min);
}



}

#endif


