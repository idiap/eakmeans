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

#ifndef ARRUTILV2L2_H
#define ARRUTILV2L2_H

#include "arrutilv2l0.h"
#include "arrutilv2l1.h"
#include "arrutilv2copy.h"
#include "arrutilv2discrete.h"


/* functions here : 
 * - generally higher level than l1
 * - set multiple variables (distances & min for example)
 * -  may use kmeans words :
 * -- centroids C 
 * -- data D 
 * -- labels L 
 * -- counts H 
 * -- sums S 
 * -- min distance dn
 * -- second lowest distance dsn
 * */

namespace arrutilv2{	
	
template <typename TInt, typename TFloat>
void set_distances_L_dn(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s,  TInt * L, TFloat * dn, TFloat * distances, TInt & n_distance_calculations){
	set_rrl2ss(ndata, dimension, data, ncentroids, C, distances, n_distance_calculations);
	set_rargminmins(ndata, ncentroids, distances,  L, dn);
}



template <typename TInt, typename TFloat>
void set_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, const TInt * const L, TFloat * const S, TInt * const H){
	std::fill_n(S, ncentroids*dimension, 0);
	std::fill_n(H, ncentroids, 0);	
	for (TInt i = 0; i < ndata; ++ i){
		++H[L[i]];
		addto(dimension, data + i*dimension, S + L[i]*dimension);
	}
}


/* the old way I did it, I get all distances first. good for blas */
template <typename TInt, typename TFloat>
void update_L_S_H_heavy(TInt ndata, TInt dimension, const TFloat * const data,  TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const S, TInt * const H, TInt & nchanges){
	std::unique_ptr<TFloat []> distances_squared (new TFloat [ndata*ncentroids]);
	set_rrl22ss(ndata, dimension, data, ncentroids, centroids, data_l22s, centroid_l22s, distances_squared.get());
	for (TInt i = 0; i < ndata; ++ i){
		TInt old_label = L[i];
		L[i] = get_argmin(ncentroids, distances_squared.get() + i*ncentroids); 
		if (old_label != L[i]){
			++H[L[i]];
			--H[old_label];
			addto(dimension, data + i*dimension, S + L[i]*dimension);
			subtractfrom(dimension, data +i*dimension, S + old_label*dimension);
			++nchanges;
		}
	}
}

/* lighter on the memory, distances calculated point by point and then discarded */
template <typename TInt, typename TFloat>
void update_L_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const S, TInt * const H, TInt & nchanges){
	std::unique_ptr<TFloat []> distances_squared (new TFloat [ncentroids]);
	
	for (TInt i = 0; i < ndata; ++ i){
		//std::cout << dimension << " " << ndata << std::endl;
		
		//std::cout << "X " << std::flush;
		
		//if (i < 1000){
			//std::cout << data_l22s[i] << " " << std::flush;
		//}
		
		//std::cout << L[i] << " " << std::flush;
		
		//std::cout << *(data + i*dimension) << " " << std::flush;
		
		//if (i == 1000){
			//std::abort();
		//}
		
		set_rl22s(dimension, data + i*dimension, ncentroids, centroids, data_l22s[i], centroid_l22s, distances_squared.get());
		TInt old_label = L[i];
		L[i] = get_argmin(ncentroids, distances_squared.get()); 
		//std::cout << i << " " << old_label << " " << L[i] << std::endl;

		if (old_label != L[i]){
			++H[L[i]];
			--H[old_label];
			addto(dimension, data + i*dimension, S + L[i]*dimension);
			subtractfrom(dimension, data +i*dimension, S + old_label*dimension);
			++nchanges;
		}
	}
	//std::cout << "done in update_L_S_H " << std::endl;
}



/* Partial template function specialisation is not trivial. See
 *  
 * http://stackoverflow.com/questions/33710309/code-duplication-prevention-2-long-functions-differing-only-in-inner-loop/
 * and
 * http://www.gotw.ca/publications/mill17.htm
 * and
 * http://artofsoftware.org/2012/12/20/c-template-function-partial-specialization/
 * 
 * 
 * but very neat way to reduce code duplication between functions. Not used here currently (see cutouts.txt)
 * 
 *   */
 
template <typename TInt, typename TFloat, bool>
struct impl_do_S_H_update_from_one{
	inline void operator()(const TInt & dimension, const TFloat * const data_i, TFloat * const S_old, const TInt & old_label, TInt * const H){
		throw std::runtime_error("in generic impl_do_S_H_subtract, should only use partial specialisations (this should be caught at compile time...)");
	}
};





template <typename TInt, typename TFloat> 
inline void update_L_S_H_batch(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const S, TInt * const H, TInt & nchanges){
	
	TInt nfullbatches = ndata/nperbatch;
	TInt nfinalbatch = ndata - nfullbatches*nperbatch;
	std::unique_ptr<TFloat []> distances_squared (new TFloat [nperbatch*ncentroids]);
	TInt old_label = 0;
	
	
	/* TFloat foo; */
	
	for (TInt bi = 0; bi < nfullbatches; ++bi){
		set_rrl22ss(nperbatch, dimension, data + bi*dimension*nperbatch, ncentroids, centroids, data_l22s +bi*nperbatch, centroid_l22s, distances_squared.get());	
		for (TInt i = nperbatch*bi; i < nperbatch*(bi + 1); ++ i){
				
			old_label = L[i];
			
			/* How can get_argmin (d) be faster than options (a), (b), (c) ? Not sure whether to compliment or insult the compiler ! I am still going with (c) as it seems like it should be at least as fast as (d) and I don't want to be a monkey at a keyboard. */
			
			
			/* (a) */
			/* TFloat * ds2 = distances_squared.get() + (i - nperbatch*bi)*ncentroids;
			 * L[i] = 0;
			 * foo = ds2[0];
			 * for (TInt j = 1; j < ncentroids; ++j){
			 * 	if (ds2[j] < foo){
			 * 		foo = ds2[j];
			 *		L[i] = j;
			 * 	}
			 * } 
			 * */
			
			/* (b) */
			/* set_argminmin(ncentroids, distances_squared.get() + (i - nperbatch*bi)*ncentroids, L[i], foo);	*/
			
			/* (c) */
			set_argmin(ncentroids, distances_squared.get() + (i - nperbatch*bi)*ncentroids, L[i]);
			
			/* (d) */
			/* L[i] = get_argmin(ncentroids, distances_squared.get() + (i - nperbatch*bi)*ncentroids); */
			
			if (old_label != L[i]){
				subtractfrom(dimension, data +i*dimension, S + old_label*dimension);
				 --H[old_label];
				++H[L[i]];
				addto(dimension, data + i*dimension, S + L[i]*dimension);
				++nchanges;
			}
		}
	}

	set_rrl22ss(nfinalbatch, dimension, data + nfullbatches*dimension*nperbatch, ncentroids, centroids, data_l22s + nfullbatches*nperbatch, centroid_l22s, distances_squared.get());	
	for (TInt i = nperbatch*nfullbatches; i < ndata; ++ i){
		TInt old_label = L[i];
		L[i] = get_argmin(ncentroids, distances_squared.get() + (i - nperbatch*nfullbatches)*ncentroids); 
		if (old_label != L[i]){
			subtractfrom(dimension, data +i*dimension, S + old_label*dimension);
			--H[old_label];
			++H[L[i]];
			addto(dimension, data + i*dimension, S + L[i]*dimension);
			++nchanges;
		}
	}
}



//as above, but now included is dn, distance to nearest. 
template <typename TInt, typename TFloat> 
inline void update_L_dn_S_H_batch(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const dn, TFloat * const S, TInt * const H, TInt & nchanges){
	
	TInt nfullbatches = ndata/nperbatch;
	TInt nfinalbatch = ndata - nfullbatches*nperbatch;
	std::unique_ptr<TFloat []> distances_squared (new TFloat [nperbatch*ncentroids]);
	TInt old_label = 0;
	//data from the full batches
	for (TInt bi = 0; bi < nfullbatches; ++bi){
		set_rrl22ss(nperbatch, dimension, data + bi*dimension*nperbatch, ncentroids, centroids, data_l22s +bi*nperbatch, centroid_l22s, distances_squared.get());
		
		for (TInt i = nperbatch*bi; i < nperbatch*(bi + 1); ++ i){
				
			old_label = L[i];
			
			set_argminmin(ncentroids, 
			distances_squared.get() + 
			(i - nperbatch*bi)*ncentroids, 
			L[i], 
			dn[i]);
			
			dn[i] = std::sqrt(std::max(static_cast<TFloat> (0), dn[i]));
			
			if (old_label != L[i]){
				subtractfrom(dimension, data +i*dimension, S + old_label*dimension);
				 --H[old_label];
				++H[L[i]];
				addto(dimension, data + i*dimension, S + L[i]*dimension);
				++nchanges;
			}
		}
	}

	set_rrl22ss(nfinalbatch, dimension, data + nfullbatches*dimension*nperbatch, ncentroids, centroids, data_l22s + nfullbatches*nperbatch, centroid_l22s, distances_squared.get());	
	for (TInt i = nperbatch*nfullbatches; i < ndata; ++ i){
		TInt old_label = L[i];
		
				
		set_argminmin(ncentroids, distances_squared.get() + (i - nperbatch*nfullbatches)*ncentroids, L[i], dn[i]);
		
		dn[i] = std::sqrt(std::max(static_cast<TFloat> (0), dn[i]));


		if (old_label != L[i]){
			subtractfrom(dimension, data +i*dimension, S + old_label*dimension);
			--H[old_label];
			++H[L[i]];
			addto(dimension, data + i*dimension, S + L[i]*dimension);
			++nchanges;
		}
	}
}




template <typename TInt, typename TFloat> 
inline void update_L_S_H_batch_increment_only(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const S, TInt * const H, TInt & nchanges){
	
	
	
	
	
	

	TInt nfullbatches = ndata/nperbatch;
	TInt nfinalbatch = ndata - nfullbatches*nperbatch;
	std::unique_ptr<TFloat []> distances_squared (new TFloat [nperbatch*ncentroids]);
	TInt old_label = 0;
	//data from the full batches
	for (TInt bi = 0; bi < nfullbatches; ++bi){
		set_rrl22ss(nperbatch, dimension, data + bi*dimension*nperbatch, ncentroids, centroids, data_l22s +bi*nperbatch, centroid_l22s, distances_squared.get());	
		for (TInt i = nperbatch*bi; i < nperbatch*(bi + 1); ++ i){
				
			old_label = L[i];
			L[i] = get_argmin(ncentroids, distances_squared.get() + (i - nperbatch*bi)*ncentroids); 
			++H[L[i]];
			addto(dimension, data + i*dimension, S + L[i]*dimension);
			if (old_label != L[i]){	
				++nchanges;
			}
		}
	}

	//data from the tail
	set_rrl22ss(nfinalbatch, dimension, data + nfullbatches*dimension*nperbatch, ncentroids, centroids, data_l22s + nfullbatches*nperbatch, centroid_l22s, distances_squared.get());	
	for (TInt i = nperbatch*nfullbatches; i < ndata; ++ i){
		TInt old_label = L[i];
		L[i] = get_argmin(ncentroids, distances_squared.get() + (i - nperbatch*nfullbatches)*ncentroids); 
		++H[L[i]];
		addto(dimension, data + i*dimension, S + L[i]*dimension);
		if (old_label != L[i]){
			++nchanges;
		}
	}
}





//As above, but now include is dn, distance to nearest
template <typename TInt, typename TFloat> 
inline void update_L_dn_S_H_batch_increment_only(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const dn, TFloat * const S, TInt * const H, TInt & nchanges){
	

	TInt nfullbatches = ndata/nperbatch;
	TInt nfinalbatch = ndata - nfullbatches*nperbatch;
	std::unique_ptr<TFloat []> distances_squared (new TFloat [nperbatch*ncentroids]);
	TInt old_label = 0;
	//data from the full batches
	for (TInt bi = 0; bi < nfullbatches; ++bi){
		set_rrl22ss(nperbatch, dimension, data + bi*dimension*nperbatch, ncentroids, centroids, data_l22s +bi*nperbatch, centroid_l22s, distances_squared.get());	
		for (TInt i = nperbatch*bi; i < nperbatch*(bi + 1); ++ i){
				
			old_label = L[i];
			set_argminmin(ncentroids, distances_squared.get() + (i - nperbatch*bi)*ncentroids, L[i], dn[i]);
			
			dn[i] = std::sqrt(std::max(static_cast<TFloat> (0), dn[i]));

			++H[L[i]];
			addto(dimension, data + i*dimension, S + L[i]*dimension);
			if (old_label != L[i]){	
				++nchanges;
			}
		}
	}

	//data from the tail
	set_rrl22ss(nfinalbatch, dimension, data + nfullbatches*dimension*nperbatch, ncentroids, centroids, data_l22s + nfullbatches*nperbatch, centroid_l22s, distances_squared.get());	
	for (TInt i = nperbatch*nfullbatches; i < ndata; ++ i){
		TInt old_label = L[i];
		//L[i] = get_argmin(ncentroids, distances_squared.get() + (i - nperbatch*nfullbatches)*ncentroids); 
		set_argminmin(ncentroids, distances_squared.get() + (i - nperbatch*nfullbatches)*ncentroids, L[i], dn[i]);
		dn[i] = std::sqrt(std::max(static_cast<TFloat> (0), dn[i]));

		++H[L[i]];
		addto(dimension, data + i*dimension, S + L[i]*dimension);
		if (old_label != L[i]){
			++nchanges;
		}
	}
}






template <typename TInt, typename TFloat>
TInt get_nchanges_update_L_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids , const TFloat * const data_l22s, const TFloat * const centroid_l22s,  TInt * const L, TFloat * const S, TInt * const H){
	TInt nchanges = 0;
	update_L_S_H(ndata, dimension, data, ncentroids, data_l22s, centroids, centroid_l22s, L, S, H, nchanges);
	return nchanges;
}

template <typename TInt, typename TFloat>
void set_L_S_H_from_distances(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const distances_rank_reserving_proxy, TInt * const L, TFloat * const S, TInt * const H){

	std::fill_n(S, ncentroids*dimension, 0);
	std::fill_n(H, ncentroids, 0);

	for (TInt i = 0; i < ndata; ++ i){
		set_argmin(ncentroids, distances_rank_reserving_proxy + i*ncentroids, L[i]);
		++H[L[i]];
		addto(dimension, data + i*dimension, S + L[i]*dimension);
	}
}
	
template <typename TInt, typename TFloat>
void set_L_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const S, TInt * const H){
	std::fill_n(S, ncentroids*dimension, 0);
	std::fill_n(H, ncentroids, 0);	
	std::unique_ptr<TFloat []> distances_squared (new TFloat [ncentroids]);
	for (TInt i = 0; i < ndata; ++ i){
		set_rl22s(dimension, data + i*dimension, ncentroids, centroids, data_l22s[i], centroid_l22s, distances_squared.get());
		L[i] = get_argmin(ncentroids, distances_squared.get()); 
		++H[L[i]];
		addto(dimension, data + i*dimension, S + L[i]*dimension);
	}
}

template <typename TInt, typename TFloat>
void set_L_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const S, TInt * const H, TInt & nchanges){
	set_L_S_H(ndata, dimension, data, ncentroids, centroids, data_l22s, centroid_l22s, L, S, H);
	nchanges += ndata;
}






template <typename TInt, typename TFloat>
/* could also be called set_rargminmin2s (here in case I search for this) ? */
void set_L2_dn(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dn, TFloat * const dsn, TInt & ndcalcs){
	std::unique_ptr<TFloat []> distances (new TFloat [ncentroids]);
	if (ncentroids < 2){
		throw std::logic_error("attempt to set_ min and argmin and min2 in set_L2_dn, but ncentroids is 0 or 1.");
	}
	for (TInt i = 0; i < ndata; ++i){
		set_rl2s(dimension, data + i*dimension, ncentroids, C, data_l22s[i], C_l22s, distances.get());
		set_argminmin2nocheck(ncentroids, distances.get(), L[i], dn[i], dsn[i]);
	}
	ndcalcs += ndata*ncentroids;
}



template <typename TInt, typename TFloat>
void update_u_deltaC_from_C_C_hist(TInt ncentroids, TInt dimension, const TFloat * const C, TInt nrounds, const TFloat * const C_hist, TFloat * const  u_deltaC, TInt & ndcalcs){
	for (TInt ri = 0; ri < nrounds; ++ri){
		set_rl2sconst(ncentroids, dimension, C, C_hist + ri*ncentroids*dimension, u_deltaC + ri*ncentroids, ndcalcs);
	}	
}



template <typename TInt, typename TFloat>
void update_C_C_l22s_from_SH(TInt ncentroids, TInt dimension, const TFloat * const S, const TInt * const H, TFloat * const C, TFloat * const C_l22s){	


	for (TInt ci = 0; ci < ncentroids; ++ci){
		if (H[ci] != 0){
			std::memcpy(C + ci*dimension, S + ci*dimension, dimension*sizeof(TFloat));
			TFloat countsi_inv = static_cast<TFloat>(1.) / static_cast<TFloat>(H[ci]);
			scale(dimension, countsi_inv, C + ci*dimension);
			set_l22(dimension, C + ci*dimension, C_l22s[ci]);
		}
		
		else {
			//let sleeping dragons lie.
		}
	}	
	//ndcalcs += ncentroids;
}

template <typename TInt, typename TFloat>
void move_then_update_C_C_l22s_from_SH(TInt ncentroids, TInt dimension, const TFloat * const S, const TInt * const H, TFloat * const C, TFloat * const C_l22s, TFloat * const C_hist, TFloat * const C_l22s_hist){
	std::memcpy(C_hist, C, ncentroids*dimension*sizeof(TFloat));
	std::memcpy(C_l22s_hist, C_l22s, ncentroids*sizeof(TFloat));
	update_C_C_l22s_from_SH(ncentroids, dimension, S, H, C, C_l22s);
}




template <typename TInt, typename TFloat>
void update_C_C_l22s_delta_C_from_SH(TInt ncentroids, TInt dimension, const TFloat * const S, const TInt * const H, TFloat * const C, TFloat * const C_l22s, TFloat * const delta_C, TInt & ndcalcs){	

	auto oldcentroids = copy_ptrarr_to_uptrarr(ncentroids*dimension, C);
	//std::fill_n(C, ncentroids*dimension, 0);
	for (TInt ci = 0; ci < ncentroids; ++ci){
		if (H[ci] != 0){
			std::memcpy(C + ci*dimension, S + ci*dimension, dimension*sizeof(TFloat));
			TFloat countsi_inv = static_cast<TFloat>(1.) / static_cast<TFloat>(H[ci]);
			scale(dimension, countsi_inv, C + ci*dimension);
		}
		
		else {
			//let setting suns set
		}
	}
	
	set_rl22s(ncentroids, dimension, C, C_l22s);
	subtractfrom(dimension*ncentroids, C, oldcentroids.get());
	set_rl2s(ncentroids, dimension, oldcentroids.get(), delta_C);
	ndcalcs += ncentroids;
}



//template <typename TInt, typename TFloat>
//void set_delta_C(const TFloat * const C_A)
	//subtractfrom(dimension*ncentroids, C, oldcentroids.get());
	//set_rl2s(ncentroids, dimension, oldcentroids.get(), delta_C);



/* see comment for set_l2gramm_partial : array pointers passed in are not dep on r0. */
template <typename TInt, typename TFloat>
void update_CC_halfminCC_partial(TInt r0, TInt r1, TInt ncentroids, TInt dimension, const TFloat * const C, const TFloat * const C_l22s,  TFloat * const CC, TFloat * const halfminCC, TInt & ndcalcs){
	set_l2gramm_partial(r0, r1, ncentroids, dimension, C, C_l22s, CC, ndcalcs);
	std::unique_ptr<TInt []> exclusions (new TInt [r1 - r0]);
	std::iota(exclusions.get(), exclusions.get() + (r1 - r0), r0);
	set_rminsexclusion(r1 - r0, ncentroids, CC + r0*ncentroids, exclusions.get(), halfminCC + r0);
	scale(r1 - r0, 0.5, halfminCC + r0);
}


/* where A and B are Nxdimension arrays, set l2s to be the distance between all rows where, def otherwise */
template <typename TFloat, typename TInt, typename Container>
void set_rl22s_where(TInt N, TInt dimension, const TFloat * const A, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, const Container & where, TFloat def, TFloat * const l22s, TInt & ndcalcs){
	std::fill_n(l22s, N, def);
	for (auto & i: where){	
		set_l22(dimension, A + i*dimension, B + i*dimension, A_l22s[i], B_l22s[i], l22s[i]);
	}
	ndcalcs += where.size();	
}

/* where A and B are Nxdimension arrays, set l2s to be the distance between all rows where, def otherwise */
template <typename TFloat, typename TInt, typename Container>
void set_rl2s_where(TInt N, TInt dimension, const TFloat * const A, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, const Container & where, TFloat def, TFloat * const l2s, TInt & ndcalcs){
	set_rl22s_where(N, dimension, A, B, A_l22s, B_l22s, where, def, l2s, ndcalcs);
	for (auto & i: where){	
		l2s[i] = std::sqrt(std::max(static_cast<TFloat>(0), l2s[i]));
	}
}

/* where A and B are Nxdimension arrays, set l2s to be the distance between all rows where, def otherwise */
template <typename TFloat, typename TInt, typename Container>
void set_rl2s_where_else_increment(TInt N, TInt dimension, const TFloat * const A, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, const Container & where, const TFloat * const incrementarr, TFloat * const l2s, TInt & ndcalcs){
	addto(N, incrementarr, l2s);
	for (auto & i: where){	
		set_l22(dimension, A + i*dimension, B + i*dimension, A_l22s[i], B_l22s[i], l2s[i]);
	}
	ndcalcs += where.size();
}




template <typename TInt, typename TFloat>
std::unique_ptr<TFloat []> get_S_from_data_and_L(TInt ndata, TInt dimension, const TFloat * const data, TInt range, const TInt * const L){
	std::unique_ptr<TFloat []> S (new TFloat[range*dimension] );
	std::fill_n(S.get(), range*dimension, 0);
	for (TInt i = 0; i < ndata; ++i){
		addto(dimension, data + i*dimension, S.get() + dimension*L[i]);
	}
	return S;
}

template <typename TInt, typename TFloat>
void set_rl22s(TInt nrows, TInt ncols, const TFloat * const A, const TFloat * const B , TFloat * const l22s, TInt & ndcals){
	std::unique_ptr<TFloat []> temp (new TFloat [nrows*ncols]);
	std::memcpy(temp.get(), A, sizeof(TFloat)*nrows*ncols);
	subtractfrom(nrows*ncols, B, temp.get());
	set_rl22s(nrows, ncols, temp.get(), l22s);
	ndcals += nrows;
}

template <typename TInt, typename TFloat>
void set_rl2s(TInt nrows, TInt ncols, const TFloat * const A, const TFloat * const B , TFloat * const l2s, TInt & ndcals){
	set_rl22s(nrows, ncols, A, B , l2s, ndcals);
	for (TInt i = 0; i < nrows; ++i){
		l2s[i] = std::sqrt(std::max(static_cast<TFloat>(0), l2s[i]));
	}
}

template <typename TInt, typename TFloat>
std::unique_ptr<TFloat []> get_rl2s(TInt nrows, TInt ncols, const TFloat * const A, const TFloat * const B, TInt & ndcalcs){
	std::unique_ptr<TFloat []> l2s (new TFloat [nrows]);
	set_rl2s(nrows, ncols, A, B, l2s.get(), ndcalcs);
	return l2s;
}


/*for each row of lowers, find for each partition the lowest value (excluding value in L)
 * TODO : make this function more efficient
 * */
template <typename TInt, typename TFloat>
void set_minlowers_excl(TInt nrows, TInt ncols, const TFloat * const lowers, const TInt * const L, TInt npartitions, TFloat * const minlowers){
	
	if (2*npartitions > ncols){
		throw std::runtime_error("2*npartitions should be greater than ncols (ncentroids?) in set_minlowers_excl");
	}
	
	TInt p_start;
	TInt p_end;
	for (TInt i = 0; i < nrows; ++i){
		for (TInt p = 0; p < npartitions; ++ p){
			minlowers[npartitions*i + p] = std::numeric_limits<TFloat>::max();
			p_start = (p*ncols)/npartitions;
			p_end = ((p+1)*ncols)/npartitions;
			for (TInt ci = p_start; ci < p_end; ++ ci){
				if (ci != L[i]){
					if (lowers[i*ncols + ci] < minlowers[npartitions*i + p]){
						minlowers[npartitions*i + p] = lowers[i*ncols + ci];
					}
				}
			}
		}
	}
}


template <typename TInt, typename TFloat>
void set_halfpartialrowsorted(TInt nrows, TInt ncols, const TFloat * const A, TInt nsorted, TFloat * const partialhalfordered, TInt * const indices_partialhalfordered){
	
	if (nsorted >= ncols){
		throw std::runtime_error("Attempt to partially sort beyond end in set_halfpartialsorted (or nsorted = ncols in which case partial sort does not make sense)");
	}
	std::vector<std::pair<TFloat, TInt>> pairs(ncols);
	
	for (TInt r = 0; r < nrows; ++ r){
		for (TInt c = 0; c < ncols; ++c){
			pairs[c].first = 0.5*A[r*ncols + c];
			pairs[c].second = c;
		}
		std::partial_sort(pairs.begin(), pairs.begin() + nsorted, pairs.end());
		for (TInt c = 0; c < nsorted; ++c){
			partialhalfordered[r*nsorted + c] = pairs[c].first;
			indices_partialhalfordered[r*nsorted + c] = pairs[c].second;
		}		
	}
}


template <typename TInt, typename TFloat>
TInt get_insertion_index(TFloat vin, TInt n_sortedv, const TFloat * const sortedv){
	TInt insertion_index = 0;
	while (insertion_index < n_sortedv){
		if (vin > sortedv[insertion_index]){
			++ insertion_index;
		}
		else{
			break;
		}
	}
	return insertion_index;
}

template <typename TInt, typename TFloat>
TInt get_insertion_index_binarial(TFloat vin, TInt n_sortedv, const TFloat * const sortedv){

	if (vin <= sortedv[0]){
		return static_cast<TInt>(0);
	}
	else if (vin > sortedv[n_sortedv - 1]){
		return  n_sortedv;
	}
	
	else{
		TInt ui = n_sortedv - 1;
		TInt li = 0;
		TInt mi;
		
		while (ui != li){
			mi = (ui + li)/2;
			if (vin > sortedv[mi]){
				li = mi +1;
			}
			else{
				ui = mi;
			}
		}
		return ui;
	}
}


/* using  indices superwhere[i0] -> superwhere[i1], set distances[0] -> distances[0 + i1 - i0] */
template <typename TInt, typename TFloat>
void set_euclidean_distances_at(TInt i0, TInt i1, const TInt * const superwhere, TInt dimension, const TFloat * const v, const TFloat * const C, TFloat * const distances, TFloat v_l22, const TFloat * const C_l22s, TInt & ndcalcs){
	
	TInt where;
	for (TInt i = 0; i < i1 - i0; ++i){
		where = superwhere[i0 + i];
		set_l2(dimension, v,  C + where*dimension, v_l22, C_l22s[where], distances[i]);		
	}
	ndcalcs += i1 - i0;
	
} 					


/* similar to set_euclidean_distances_at  */
template <typename TInt, typename TFloat>
void set_l2s_at(TInt nwhere, const TInt * const allwhere, TInt dimension, const TFloat * const v, const TFloat * const C, TFloat v_l22, const TFloat * const C_l22s, TFloat * const distances, TInt & ndcalcs){
	
	
	TInt where;
	for (TInt i = 0; i < nwhere; ++i){
		where = allwhere[i];
		set_l2(dimension, v,  C + where*dimension, v_l22, C_l22s[where], distances[i]);		
	}
	ndcalcs += nwhere;	
}






/* set glowers, dn, L, group, */
template <typename TInt, typename TFloat>
void set_L_group_dn_glowers(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const C_l22s,  TInt ngroups, const TInt * const groupparts, const TInt * const groupsizes, TInt * const  L, TInt * const group, TFloat * const  dn, TFloat * const  glowers, TInt & ndcalcs){
	
	
	
	
	std::unique_ptr<TFloat []> uptr_distances (new TFloat [ncentroids]);
	TFloat * const distances = uptr_distances.get();
	
	TInt label_group_nearest; 
	TFloat distance_group_nearest;
	TFloat distance_group_second_nearest;
	
	for (TInt i = 0; i < ndata; ++i){
		dn[i] = std::numeric_limits<TFloat>::max();
		set_rl2s(dimension, data + i*dimension, ncentroids, C, data_l22s[i], C_l22s, distances, ndcalcs);
		group[i] = 0;
		for (TInt gi = 0; gi < ngroups; ++gi){
			set_argminmin2(groupsizes[gi], distances + groupparts[gi], label_group_nearest, distance_group_nearest, distance_group_second_nearest);
			
			if (distance_group_nearest < dn[i]){
				//if gi != group[i], this sets the lower bound in group[i]. if gi == group[i], (then gi = 0) : no effect.
				glowers[i*ngroups + group[i]] = dn[i];
				L[i] = groupparts[gi] + label_group_nearest;
				dn[i] = distance_group_nearest;
				glowers[i*ngroups + gi] = distance_group_second_nearest;
				group[i] = gi;
			}
			else{
				glowers[i*ngroups + gi] = distance_group_nearest;
			}
		}
	}
	

	
}

template <typename TInt, typename TFloat>
void update_delta_G(TInt ngroups, const TFloat * const delta_C, const TInt * const groupparts, TFloat * const delta_G){
	for (TInt gi = 0; gi < ngroups; ++gi){
		delta_G[gi] = delta_C[groupparts[gi]];
		for (TInt ci = groupparts[gi] + 1; ci <  groupparts[gi + 1]; ++ci){
			if (delta_C[ci] > delta_G[gi]){
				delta_G[gi] = delta_C[ci];
			}
		}
	}
}

	

//does this function already exist under another name ? 
template <typename TInt, typename TFloat>
void set_rl2s_at_L(TInt dimension, TInt ndata, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TInt * const L, TFloat * const dn){  
	for (TInt i = 0; i < ndata; ++i){
		set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], dn[i]);
	}
}


template <typename TInt, typename TFloat>
void set_rl2s_at_L(TInt dimension, TInt ndata, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TInt * const L, TFloat * const dn, TInt & ndcalcs){  
	for (TInt i = 0; i < ndata; ++i){
		set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], dn[i]);
	}
	ndcalcs += ndata;
}


/* variance  of ndata points in dimension dimension*/
template <typename TInt, typename TFloat>
TFloat get_variance(TInt dimension, TInt ndata, const TFloat * const data, const TFloat * const data_l22s){
	TFloat sum_l22s = 0;
	for (TInt i = 0; i < ndata; ++ i){
		sum_l22s += data_l22s[i];
	}
	
	std::unique_ptr<TFloat []> sum_data (new TFloat [dimension]);
	set_sums(ndata, dimension, data, sum_data.get(), false);
	
	TFloat the_smaller_term = 0;
	for (TInt di = 0; di < dimension; ++di){
		the_smaller_term += sum_data[di]*(sum_data[di]/static_cast<TFloat>(ndata));
	}
	
	return sum_l22s - the_smaller_term;
}

template <typename TInt, typename TFloat>
inline void set_rl2s_argminmin(const TInt & dimension, const TFloat * const v, const TInt & nrowsC, const TFloat * const C, const TFloat & v_l22, const TFloat * const C_l22s, TFloat * const l2s, TInt & argmin, TFloat & min){
	set_rl2s(dimension, v, nrowsC, C, v_l22, C_l22s, l2s);
	set_argminmin(nrowsC, l2s, argmin, min);	
}

template <typename TInt, typename TFloat>
void set_rrl2ss_argminmins(TInt nrowsA, TInt dimension, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l2ss, TInt * const argmins, TFloat * const mins){
	
		
	set_rrl22ss(nrowsA, dimension, A, nrowsB, B, A_l22s, B_l22s, l2ss);
	for (TInt i = 0; i < nrowsA; ++i){
		for (TInt j = 0; j < nrowsB; ++j){
			l2ss[nrowsB*i + j] = std::sqrt(std::max(static_cast<TFloat>(0), l2ss[nrowsB*i + j]));
		}
		set_argminmin(nrowsB, l2ss + nrowsB*i, argmins[i], mins[i]);
	}		
}

template <typename TInt, typename TFloat>
void set_rrl2ss_argminmins(TInt nrowsA, TInt dimension, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l22ss, TInt * const argmins, TFloat * const mins, TInt & ndcalcs){
	set_rrl2ss_argminmins(nrowsA, dimension, A, nrowsB, B, A_l22s, B_l22s, l22ss, argmins,  mins);
	ndcalcs += nrowsA*nrowsB;
}
 

template <typename TInt, typename TFloat>
inline void  update_delta_G_from_delta_C(const TInt & ncentroids, const TFloat * const delta_C, const TInt &  ngroups, const TInt * const groupparts, TFloat * const delta_G){
	for (TInt gi = 0; gi < ngroups; ++gi){
		delta_G[gi] = delta_C[groupparts[gi]];
		for (TInt ci = groupparts[gi] + 1; ci < groupparts[gi+1]; ++ci){
			if (delta_C[ci] > delta_G[gi]){
				delta_G[gi] = delta_C[ci];
			}
		}
	}
}

template <typename TInt, typename TFloat>
void update_u_delta_G_global_from_u_delta_C(TInt nrounds, TInt ncentroids, const TFloat * const u_delta_C, TInt ngroups, const TInt * const groupparts, TFloat * const u_delta_G, TFloat * const u_delta_global){
	for (TInt r = 0; r< nrounds; ++r){
		update_delta_G_from_delta_C(ncentroids, u_delta_C + r*ncentroids, ngroups, groupparts, u_delta_G + r*ngroups);
		u_delta_global[r] = *std::max_element(u_delta_G + r*ngroups, u_delta_G + (r + 1)*ngroups );
	}
}

/* déja vu ? */
template <typename TInt, typename TFloat>
void update_cumabs_max(TInt ncentroids, const TFloat * const delta_C, TInt nrounds, TFloat * const cumabs, TFloat * const max_deltaC_since){
	for (TInt r = 0; r < nrounds; ++r){
		addto(ncentroids, delta_C, cumabs + r*ncentroids);
		/* My inner speed freak says shouldn't get max like this (why monitor argmax if you don't need to?) */
		max_deltaC_since[r] = * std::max_element(cumabs + r*ncentroids, cumabs + (r + 1)*ncentroids); 
	}
}


//Sets L, dn and lowers and increments S and H.
template <typename TInt, typename TFloat> 
inline void set_L_lowers_dn_and_increment_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt * const L, TFloat * const lowers, TFloat * const dn, TFloat * const S, TInt * const H, TInt & nchanges){

	//set L, dn and lowers
	set_rrl2ss_argminmins(ndata, dimension, data, ncentroids, centroids, data_l22s, centroid_l22s, lowers, L, dn);
	
	
	//TInt nrowsA, TInt dimension, const TFloat * const A, TInt nrowsB, const TFloat * const B, const TFloat * const A_l22s, const TFloat * const B_l22s, TFloat * const l2ss, TInt * const argmins, TFloat * const mins
	
	//increment S, H.	
	for (TInt i = 0; i < ndata; ++ i){
		++H[L[i]];
		addto(dimension, data + i*dimension, S + L[i]*dimension);
	}
	nchanges += ndata;

}









}



#endif
