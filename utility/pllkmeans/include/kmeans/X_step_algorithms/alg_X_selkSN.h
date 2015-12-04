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

#ifndef ALG_X_SELKSN_H
#define ALG_X_SELKSN_H

#include "arrutilv2l1.h"
#include "arrutilv2l2.h"
#include "arrutilv2l3.h"

namespace kmeans{

template <typename TInt, typename TFloat>
/* technical note : round is not used in this function, it is a vestige of video making */
void update_L_lowers_upbs_S_H_3v0(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt &ndcalcs, 
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

}


#endif
