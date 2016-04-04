#ifndef ARRUTILV2L0BLASLESS_H
#define ARRUTILV2L0BLASLESS_H

namespace arrutilv2{

inline void proxy_openblas_set_num_threads(int nthreads){

}


template <typename TInt, typename TFloat>
void set_rargmaxabs(TInt nrows, TInt ncols, const TFloat * const A,  TInt * argmaxabss){
	TFloat absval;
	TFloat vmax;
	for (TInt r = 0; r < nrows; ++r){
		argmaxabss[r] = 0;
		vmax = A[r*ncols];
		for (TInt c = 1; c < ncols; ++c){
			absval = std::abs(A[r*ncols + c]);
			if (absval > vmax){
				argmaxabss[r] = c;
				vmax = absval;
			}
		}
	}
}

template <typename TInt, typename TFloat, typename TScaleFloatType>
void scale(TInt N, TScaleFloatType factor, TFloat * const toscale){
	TFloat factor_ct = static_cast<TFloat> (factor);
	for (TInt i = 0; i < N; ++i){
		toscale[i] *= factor_ct;
	}
}



template <typename TInt, typename TFloat>
void rank1rowupdate(TInt ncols,  const TFloat * const row, TFloat scale, TInt nrows, TFloat * const toupdate){
	std::unique_ptr<TFloat []> scaledrow (new TFloat [ncols]);
	for (TInt c = 0; c < ncols; ++c){
		scaledrow[c] = row[c]*scale;
	}
	for (TInt r = 0; r < nrows; ++r){
		for (TInt c = 0; c < ncols; ++c){
			toupdate[r*ncols + c] += scaledrow[c];
		}
	}
}		

template <typename TInt, typename TFloat>
inline void set_l22(const TInt & ndata, const TFloat * const a, TFloat & l22){
	l22 = 0;
	for (TInt i = 0; i < ndata; ++i){
		l22 += a[i]*a[i];
	}
}

template <typename TInt, typename TFloat>
inline void set_l22(const TInt & dimension, const TFloat * const a, const TFloat * const b, const TFloat & a_l22, const TFloat & b_l22, TFloat & l22){
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
	sum = 0;
	for (TInt i = 0; i < ndata; ++i){
		sum += a[i];
	}
}

template <typename TInt, typename TFloat>
void set_l22s(TInt nrows, TInt ncols, const TFloat * const A, TFloat * const l22s, bool byrow){
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
	for (TInt r = 0; r < nrows; ++r){
		l22s[r] = 0;
		for (TInt c = 0; c < ncols; ++c){
			l22s[r] += B[r*ncols + c]*v[c];
		}
		l22s[r] *= -2;
		l22s[r] += v_l22s;
		l22s[r] += B_l22s[r];
	}
}


/* A : nrowsA x ncols
* B : nrowsB x ncols
* C[i,j] : |A[i] - B[j]|_2 for 0 <= i < nrowsA   0 <= j < nrowsB 
*  */

template <typename TInt, typename TFloat>
void set_rrl22ss(TInt nrowsA, TInt ncols, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l22ss){
	for (TInt r = 0; r < nrowsA; ++r){
		set_rl22s(ncols, A + r*ncols, nrowsB, B, A_l22s[r], B_l22s, l22ss + r*nrowsB);
	}
}


template <typename TInt, typename TFloat>
void subtractfrom(TInt N, const TFloat * const tosubtract, TFloat * const from){
	for (TInt i = 0; i < N; ++i){
		from[i] -= tosubtract[i];
	}
}	

template <typename TInt, typename TFloat>
void addto(TInt N, const TFloat * const toadd, TFloat * const to){
	for (TInt i = 0; i < N; ++i){
		to[i] += toadd[i];
	}
}





}

#endif




