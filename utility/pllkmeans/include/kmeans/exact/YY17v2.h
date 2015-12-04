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

#ifndef PLL_ELKANKMEANS_17V2_H
#define PLL_ELKANKMEANS_17V2_H

#include "baseYYSMN.h" //for initial clustering of the centroids


#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{
	

template <typename TInt, typename TFloat>
/* the same as a17v2.h : no delta_C test as per yy, so if group bound fails all distances calculated. 
 * a few percent slower, strange as most pxx are faster than axx.
 * */
void update_L_glowers_upb_S_H_17v2(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt & ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const  C_l22s, const TFloat * const delta_C, const TFloat * const  delta_G,  TInt * const  L, TInt ngroups, const TInt * const groupparts, const TInt * const  groupsizes, TInt * const  group, TFloat * const glowers, TFloat * const upb, const TInt & round){
	

	std::unique_ptr<TFloat []> distances( new TFloat [ncentroids] );	
	TInt group_nearest_index;
	TFloat group_nearest;
	TFloat group_second_nearest;

	for (TInt i = 0; i < ndata; ++i){
		arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upb[i] , ndcalcs );
		TInt label_before = L[i];	
			
		for (TInt gi = 0; gi < ngroups; ++gi){
			glowers[i*ngroups + gi] -= delta_G[gi];
		}
		
		for (TInt gi = 0; gi < ngroups; ++gi){
			if (glowers[i*ngroups + gi] < upb[i]){
				arrutilv2::set_rl2s(dimension, data + i*dimension, groupsizes[gi], C + groupparts[gi]*dimension, data_l22s[i], C_l22s + groupparts[gi], distances.get(), ndcalcs);

				if (gi != group[i]){
					arrutilv2::set_argminmin2(groupsizes[gi], distances.get(), group_nearest_index, group_nearest, group_second_nearest);
					group_nearest_index += groupparts[gi];								
					if (group_nearest < upb[i]){
						if (gi < group[i]){
							glowers[i*ngroups + group[i]] = std::min(upb[i], glowers[i*ngroups + group[i]]);
						}
						else{
							glowers[i*ngroups + group[i]] = upb[i];
						}
						glowers[i*ngroups + gi] = group_second_nearest;
						L[i] = group_nearest_index;
						group[i] = gi;
						upb[i] = group_nearest;
					}

					else{
						glowers[i*ngroups + gi] = group_nearest;
					}
				}

				else{
					arrutilv2::set_argminmin2(groupsizes[gi], distances.get(), L[i], upb[i], glowers[i*ngroups + gi]);
					L[i] += groupparts[gi];
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





template <typename TInt, typename TFloat>
class P17V2 : public kmeans::BaseYYSMN<TInt, TFloat>{
	

	protected:	
		
		std::function<void(TInt)> update_L_glowers_upb_S_H_17v2_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				this->pll_principal_X(update_L_glowers_upb_S_H_17v2<TInt, TFloat>, ti, 
				this->get_delta_C(), this->get_delta_G(), this->get_L() + x0, this->get_ngroups(), this->get_groupparts(), this->get_groupsizes(), this->get_group() + x0, this->get_glowers() + x0*this->get_ngroups(), this->get_upb() + x0, this->round);
			};
		}
		
	public:
		typedef kmeans::BaseYYSMN<TInt, TFloat> BC;
		template<typename... Args>
		P17V2(Args&&... args): BC(std::forward<Args>(args)...)
		
		{
			this->setalgname("p17v2");
		}
		
		virtual ~P17V2(){}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to 17v2..\n\n";
		}

		virtual void set_initialisation_tasks(){
			this->yinyang_initialisation_tasks();
		}
		
		virtual void set_C_tasks(){
			
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX),

				arrutilv2::update_delta_G_from_delta_C_ati(this->getncentroids(), this->get_delta_C(), this->get_ngroups(), this->get_groupparts(), this->get_groupsizes(), this->get_delta_G())
			};
		}
		
		virtual void set_X_tasks(){
			
			this->X_tasks = {
				this->update_L_glowers_upb_S_H_17v2_ati()
			};
		}
};

}

#endif

