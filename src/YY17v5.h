/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_ELKANKMEANS_17V5_H
#define PLL_ELKANKMEANS_17V5_H

#include "baseYYSMN.h"


#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{
	

template <typename TInt, typename TFloat>
/* exactly a17v5 (yinyang with lowerbound set to distance immediately)
 * */

void update_L_glowers_upb_S_H_17v5(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt & ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const  C_l22s, const TFloat * const delta_C, const TFloat * const  delta_G,  TInt * const  L, TInt ngroups, const TInt * const groupparts, const TInt * const  groupsizes, TInt * const  group, TFloat * const glowers, TFloat * const upb, const TInt & round){
	
	/* the local workforce */
	std::vector<TFloat> glowers_previous(ngroups);	
	TFloat adist;
	TFloat group_nearest;
	TFloat group_second_nearest;
	TFloat group_nearest_index;
	std::unique_ptr<TFloat []> distances (new TFloat [ncentroids]);	
	TFloat globallower;
	
	for (TInt i = 0; i < ndata; ++i){

		TInt label_before = L[i];
		globallower = std::numeric_limits<TFloat>::max();			
		for (TInt gi = 0; gi < ngroups; ++gi){
			glowers_previous[gi] = glowers[i*ngroups + gi];
			glowers[i*ngroups + gi] -= delta_G[gi];
			globallower =std::min(globallower, glowers[i*ngroups + gi]);
		}
		
		arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]],  upb[i], ndcalcs);
		if (globallower < upb[i]){
			for (TInt gi = 0; gi < ngroups; ++gi){
				if (glowers[i*ngroups + gi] < upb[i]){
					
					group_nearest =  std::numeric_limits<TFloat>::max();	
					group_second_nearest = std::numeric_limits<TFloat>::max();
					group_nearest_index = 0;
	
					for (TInt ci = groupparts[gi]; ci < groupparts[gi + 1]; ++ci){
						//glowers_previous[gi] - delta_C[ci] is a lower bound on ci iff ci != label_before.
						if (label_before == ci || glowers_previous[gi] - delta_C[ci] < group_second_nearest){
							arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], adist, ndcalcs);
							
							if (adist < group_nearest){
								group_second_nearest = group_nearest;
								group_nearest = adist;
								group_nearest_index = ci;
							}
							
							else if (adist < group_second_nearest){
								group_second_nearest = adist;
							}
						}
					}
						
					
					if (group[i] != gi){
						if (group_nearest < upb[i]){
							if (gi < group[i]){
								glowers[i*ngroups + group[i]] = std::min(upb[i], glowers[i*ngroups + group[i]]);
							}
							else{// group[i] < gi
								glowers[i*ngroups + group[i]] = upb[i];
							}
							glowers[i*ngroups + gi] = group_second_nearest;
							group[i] = gi;
							L[i] = group_nearest_index;
							upb[i] = group_nearest;
						}
						
						else{
							glowers[i*ngroups + gi] = group_nearest;
						}
					}
					
					else{
						glowers[i*ngroups + gi] = group_second_nearest;
						upb[i] = group_nearest;
						L[i] = group_nearest_index;
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
}





template <typename TInt, typename TFloat>
class P17V5 : public kmeans::BaseYYSMN<TInt, TFloat>{
	

	protected:	
		
		std::function<void(TInt)> update_L_glowers_upb_S_H_17v5_ati(){

			return [this](TInt ti){
		
				TInt x0 = (ti*this->getndata())/this->getnthreads();			
				this->pll_principal_X(update_L_glowers_upb_S_H_17v5<TInt, TFloat>, ti, this->get_delta_C(), this->get_delta_G(), this->get_L() + x0, this->get_ngroups(), this->get_groupparts(), this->get_groupsizes(), this->get_group() + x0, this->get_glowers() + x0*this->get_ngroups(), this->get_upb() + x0, this->round);
	
			};
		}
		
	public:
		typedef kmeans::BaseYYSMN<TInt, TFloat> BC;
		template<typename... Args>
		P17V5(Args&&... args): BC(std::forward<Args>(args)...)
		
		{
			this->setalgname("p17v5");
		}
		
		virtual ~P17V5(){}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to 17v5..\n\n";
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
				this->update_L_glowers_upb_S_H_17v5_ati()
			};
		}
};

}

#endif
