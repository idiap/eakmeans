/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef ALG_X_SELKSN_H
#define ALG_X_SELKSN_H


#include "arrutilv2l1.h"
#include "arrutilv2l2.h"
#include "arrutilv2l3.h"
#include "sparsedatasets.h"
#include "sparseutil.h"

namespace kmeans{

template <typename TInt, typename TFloat>
/* technical note : round is not used in this function, it is a vestige of video making */
void update_L_lowers_upper_S_H_3v0(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const delta_C, TInt * const L, TFloat * const lowers, TFloat * const upbs, const TInt & round){
	
	nchanges = 0;
	ndcalcs = 0;
	arrutilv2::rank1rowupdate(ncentroids, delta_C, static_cast<TFloat>(-1.), ndata, lowers);
	
	for (TInt i = 0; i < ndata; ++i){

		//update upper bound
		upbs[i] += delta_C[L[i]];
		TInt label_before = L[i];

		//the upperbound is initially not tight, requiring the behaviour in first while loop
		TInt ci = 0;
		while (ci < ncentroids){
			
			//local test with loose upperbound
			if  ((L[i] != ci) && (upbs[i] > lowers[i*ncentroids + ci])){
				
				//update loose upperbound to true distance
				arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upbs[i], ndcalcs);

				lowers[i*ncentroids + L[i]] = upbs[i];
				//same test with tight upperbound
				if  ((upbs[i] > lowers[i*ncentroids + ci])){
					arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lowers[i*ncentroids + ci], ndcalcs);
					if (upbs[i] > lowers[i*ncentroids + ci]){
						upbs[i] = lowers[i*ncentroids + ci];
						L[i] = ci;
					}
					
				}
				++ci;
				break;
			}
			++ci;
		}
		
		while (ci < ncentroids){
			//local test with tight upperbound (made tight by some preceding centroid test)
			if ((upbs[i] > lowers[i*ncentroids + ci]) ){
				arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lowers[i*ncentroids + ci], ndcalcs);
				if (upbs[i] > lowers[i*ncentroids + ci]){
					upbs[i] = lowers[i*ncentroids + ci];
					L[i] = ci;
				}
			}
			++ci;
		}
		
		
		if (L[i] != label_before){
			++nchanges;
			++H[L[i]];
			--H[label_before];
			arrutilv2::addto(dimension, data + i*dimension, S + dimension*L[i]);
			arrutilv2::subtractfrom(dimension, data + i*dimension, S + dimension*label_before);
		}
	}
}




//sparse version of 3v0 update (ideas for merging?)
template <typename TInt, typename TFloat>
void sparse_update_L_lowers_upper_where_changes_3v0(TInt ncentroids, TInt dimension, TInt data0, TInt data1, const sparse::SparseData<TInt, TFloat> & data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const delta_C, std::vector<std::tuple<TInt, TInt, TInt > > & where_label_change, TInt &ndcalcs, TInt * const L, TFloat * const lowers, TFloat * const upbs){
	
	
	where_label_change.clear();
	TInt ndata = data1 - data0;	
	ndcalcs = 0;
	arrutilv2::rank1rowupdate(ncentroids, delta_C, static_cast<TFloat>(-1.), ndata, lowers);
	
	TInt datai = 0;
	for (TInt i = 0; i < ndata; ++i){
		datai = data0 + i;
		
		upbs[i] += delta_C[L[i]];
		TInt label_before = L[i];

		TInt ci = 0;
		while (ci < ncentroids){
			if  ((L[i] != ci) && (upbs[i] > lowers[i*ncentroids + ci])){
				sparse::set_l2( 
				data.starts[datai+1] - data.starts[datai], data.indices.data() +  data.starts[datai], data.values.data() +  data.starts[datai], C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upbs[i], ndcalcs);
				
				lowers[i*ncentroids + L[i]] = upbs[i];
				if  ((upbs[i] > lowers[i*ncentroids + ci])){
					sparse::set_l2(data.starts[datai+1] - data.starts[datai], data.indices.data() +  data.starts[datai], data.values.data() +  data.starts[datai], C + ci*dimension, data_l22s[i], C_l22s[ci], lowers[i*ncentroids + ci], ndcalcs);
				
					if (upbs[i] > lowers[i*ncentroids + ci]){
						upbs[i] = lowers[i*ncentroids + ci];
						L[i] = ci;
					}
					
				}
				++ci;
				break;
			}
			++ci;
		}
		
		while (ci < ncentroids){
			if ((upbs[i] > lowers[i*ncentroids + ci]) ){
				sparse::set_l2(data.starts[datai+1] - data.starts[datai], data.indices.data() +  data.starts[datai], data.values.data() +  data.starts[datai], C + ci*dimension, data_l22s[i], C_l22s[ci], lowers[i*ncentroids + ci], ndcalcs);				
				if (upbs[i] > lowers[i*ncentroids + ci]){
					upbs[i] = lowers[i*ncentroids + ci];
					L[i] = ci;
				}
			}
			++ci;
		}
		
		if (L[i] != label_before){
			where_label_change.emplace_back(datai, label_before, L[i]);
		}
	}
}

	


template <typename TInt, typename TFloat>
/* upbs[i] on entry is distance to currently assigned nearest neighbor. On exit it is distance to now exact (and thus assigned) nearest neighbor. Note : same signature as 3v0 version.  */
void update_L_lowers_upper_S_H_3v1(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const delta_C, TInt * const L, TFloat * const lowers, TFloat * const upbs, const TInt & round){
	
	
	//std::vector<TFloat> v_l2s(ncentroids);
	
	nchanges = 0;
	ndcalcs = 0;
	arrutilv2::rank1rowupdate(ncentroids, delta_C, static_cast<TFloat>(-1.), ndata, lowers);
	TInt label_before;
	for (TInt i = 0; i < ndata; ++i){
		label_before = L[i];
		arrutilv2::set_l2(dimension, data + i*dimension, C + dimension*L[i], data_l22s[i], C_l22s[L[i]], upbs[i]);
		for (TInt ci = 0; ci < ncentroids; ++ci){
			if ((upbs[i] > lowers[i*ncentroids + ci]) ){
				arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lowers[i*ncentroids + ci], ndcalcs);
				if (upbs[i] > lowers[i*ncentroids + ci]){
					upbs[i] = lowers[i*ncentroids + ci];
					L[i] = ci;
				}
			}
		}
		
		if (L[i] != label_before){
			++nchanges;
			++H[L[i]];
			--H[label_before];
			arrutilv2::addto(dimension, data + i*dimension, S + dimension*L[i]);
			arrutilv2::subtractfrom(dimension, data + i*dimension, S + dimension*label_before);
		}
	}
}





////inner between sparse and dense
//template <typename TInt, typename TFloat>
//inline TFloat get_inner( 
//TInt n_a, const TInt * const  a_indices, const TFloat *  const a_values,
//const TFloat * const b 
//){
	
//}


/* sparse version of 3v1, identical signature as 3v0 */
template <typename TInt, typename TFloat>
void sparse_update_L_lowers_upper_where_changes_3v1(TInt ncentroids, TInt dimension, TInt data0, TInt data1, const sparse::SparseData<TInt, TFloat> & data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const delta_C, std::vector<std::tuple<TInt, TInt, TInt > > & where_label_change, TInt &ndcalcs, TInt * const L, TFloat * const lowers, TFloat * const upbs){
	
	where_label_change.clear();
	TInt ndata = data1 - data0;	
	ndcalcs = 0;
	
	arrutilv2::rank1rowupdate(ncentroids, delta_C, static_cast<TFloat>(-1.), ndata, lowers);
	TInt datai = 0; 
	TInt label_before;
	for (TInt i = 0; i < ndata; ++i){
		datai = data0 + i; 
		label_before = L[i];
		//set upbs[i] to distance from data i to centroid L[i].
		sparse::set_l2(data.starts[datai+1] - data.starts[datai], data.indices.data() +  data.starts[datai], data.values.data() +  data.starts[datai], C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upbs[i], ndcalcs);
		
	
		
		for (TInt ci = 0; ci < ncentroids; ++ci){
			if ((upbs[i] > lowers[i*ncentroids + ci]) ){
				//set lowers [i, ci] to distance from data i to centroid ci.
				sparse::set_l2(data.starts[datai+1] - data.starts[datai], data.indices.data() +  data.starts[datai], data.values.data() +  data.starts[datai], C + ci*dimension, data_l22s[i], C_l22s[ci], lowers[i*ncentroids + ci], ndcalcs);
				
				
				
				
				if (upbs[i] > lowers[i*ncentroids + ci]){
					upbs[i] = lowers[i*ncentroids + ci];
					L[i] = ci;
				}
			}
		}
		if (L[i] != label_before){
			where_label_change.emplace_back(datai, label_before, L[i]);
		}
	}
}




}


#endif
