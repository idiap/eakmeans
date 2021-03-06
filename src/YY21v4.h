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

#ifndef PLL_ELKANKMEANS_21V4_H
#define PLL_ELKANKMEANS_21V4_H

#include "baseYYMNS.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{
	

/* like 21v3, but upb_at_last in NS-ised. If anything, marginally slower than 21v3.
 * */
 
/* note that 17v3 may have fewer distance calculations due to it's (not so?) silly way of getting a global lower bound. */
 
template <typename TInt, typename TTau, typename TFloat>
 void update_L_glowers_upb_S_H_21v4_not_modround(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt & ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const  C_l22s, const TFloat * const u_delta_C, const TFloat * const u_delta_G, const TFloat * const u_delta_global, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt ngroups, const TInt * const groupparts, const TInt * const  groupsizes, const TInt & modround, TInt * const  group, TInt * const L, TTau * const tau, TTau * const tau_globallowers, TFloat * const glowers_at_last, TFloat * const globallowers_at_last,  TFloat * const upb_at_last, TTau * const tau_ub){
	
	/* some worker bees */
	std::unique_ptr<TFloat []> distances ( new TFloat [ncentroids] );	
	TInt group_nearest_index;
	TFloat group_nearest;
	TFloat group_second_nearest;	
	TFloat globallower;
	TFloat upperbound;
	for (TInt i = 0; i < ndata; ++i){
		/* TODO : using SN in NS scheme is not so cool, consider NS-ising upb (make it upb_at_last) 
		 * TOCONSIDER : do we really need u_delta_C ? Would we not be better off using it only locally while updating delta_G ? */
		
		upperbound = upb_at_last[i] + u_delta_C[ncentroids*(tau_ub[i]) + L[i]];
		//upb_at_last[i] += u_delta_C[ncentroids*(modround - 1) + L[i]];
		
		globallower = globallowers_at_last[i] - u_delta_global[tau_globallowers[i]];
		
		//if (globallower < upb_at_last[i]){
		if (globallower < upperbound){
			arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upb_at_last[i], ndcalcs);
			tau_ub[i] = modround;
			if (globallower < upb_at_last[i]){
				TInt label_before = L[i];
				globallowers_at_last[i] = std::numeric_limits<TFloat>::max();
				for (TInt gi = 0; gi < ngroups; ++gi){
					if (glowers_at_last[i*ngroups + gi] - u_delta_G[ngroups*tau[i*ngroups + gi] + gi] < upb_at_last[i]){						
						tau[i*ngroups + gi] = modround;
						arrutilv2::set_rl2s(dimension, data + i*dimension, groupsizes[gi], C + groupparts[gi]*dimension, data_l22s[i], C_l22s + groupparts[gi], distances.get(), ndcalcs);
						if (gi != group[i]){
							arrutilv2::set_argminmin2(groupsizes[gi], distances.get(), group_nearest_index, group_nearest, group_second_nearest);
							group_nearest_index += groupparts[gi];
							if (group_nearest < upb_at_last[i]){
								if (gi < group[i]){
									glowers_at_last[i*ngroups + group[i]] = std::min(upb_at_last[i],
									glowers_at_last[i*ngroups + group[i]] - u_delta_G[ngroups*tau[i*ngroups + group[i]] + group[i]]);								
								}
								else{
									glowers_at_last[i*ngroups + group[i]] = upb_at_last[i];
									globallowers_at_last[i] = upb_at_last[i];
								}
								tau[i*ngroups + group[i]] = modround;
								glowers_at_last[i*ngroups + gi] = group_second_nearest;
								L[i] = group_nearest_index;
								group[i] = gi;
								upb_at_last[i] = group_nearest;
							}
							else{
								glowers_at_last[i*ngroups + gi] = group_nearest;
							}
						}					
						else{
							arrutilv2::set_argminmin2(groupsizes[gi], distances.get(), L[i], upb_at_last[i], glowers_at_last[i*ngroups + gi]);
							L[i] += groupparts[gi];
						}
					}
					globallowers_at_last[i] = std::min(globallowers_at_last[i], glowers_at_last[i*ngroups + gi] - u_delta_G[ngroups*tau[i*ngroups + gi] + gi]);
				}
				tau_globallowers[i] = modround;
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
}


template <typename TInt, typename TTau, typename TFloat>
 void update_L_glowers_upb_S_H_21v4_modround(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H, TInt & nchanges, TInt & ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const  C_l22s, const TFloat * const u_delta_C, const TFloat * const u_delta_G, const TFloat * const u_delta_global, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt ngroups, const TInt * const groupparts, const TInt * const  groupsizes, const TInt & modround, const TInt & cutperiod, TInt * const  group, TInt * const L, TTau * const tau, TTau * const tau_globallowers, TFloat * const glowers_at_last, TFloat * const globallowers_at_last,  TFloat * const upb_at_last, TTau * const tau_ub){
	
	/* some worker bees */
	std::unique_ptr<TFloat []> distances ( new TFloat [ncentroids] );	
	TInt group_nearest_index;
	TFloat group_nearest;
	TFloat group_second_nearest;	
	
	for (TInt i = 0; i < ndata; ++i){
		
	
		for (TInt gi = 0; gi < ngroups; ++gi){
			glowers_at_last[i*ngroups + gi] -= u_delta_G[ngroups*tau[i*ngroups + gi] + gi];
			
		}
		globallowers_at_last[i] -= u_delta_global[tau_globallowers[i]];
		std::fill_n(tau + i*ngroups, ngroups, static_cast<TTau>(0));
		upb_at_last[i] += u_delta_C[ncentroids*(tau_ub[i]) + L[i]];
		tau_ub[i] = 0;
		tau_globallowers[i] = 0;


		//upb_at_last[i] += u_delta_C[ncentroids*(cutperiod - 1) + L[i]];
		
		if (globallowers_at_last[i] < upb_at_last[i]){
		//if (globallowers_at_last[i] < upperbound){
			arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upb_at_last[i], ndcalcs);
			tau_ub[i] = modround;
			if (globallowers_at_last[i] < upb_at_last[i]){
				TInt label_before = L[i];
				globallowers_at_last[i] = std::numeric_limits<TFloat>::max();
				for (TInt gi = 0; gi < ngroups; ++gi){
					if (glowers_at_last[i*ngroups + gi] < upb_at_last[i]){						
						//tau[i*ngroups + gi] = modround;
						arrutilv2::set_rl2s(dimension, data + i*dimension, groupsizes[gi], C + groupparts[gi]*dimension, data_l22s[i], C_l22s + groupparts[gi], distances.get(), ndcalcs);
						if (gi != group[i]){
							arrutilv2::set_argminmin2(groupsizes[gi], distances.get(), group_nearest_index, group_nearest, group_second_nearest);
							group_nearest_index += groupparts[gi];
							if (group_nearest < upb_at_last[i]){
								if (gi < group[i]){
									glowers_at_last[i*ngroups + group[i]] = std::min(upb_at_last[i],
									glowers_at_last[i*ngroups + group[i]]);// - u_delta_G[ngroups*tau[i*ngroups + group[i]] + group[i]]);								
								}
								else{
									glowers_at_last[i*ngroups + group[i]] = upb_at_last[i];
									globallowers_at_last[i] = upb_at_last[i];
								}
								//tau[i*ngroups + group[i]] = modround;
								glowers_at_last[i*ngroups + gi] = group_second_nearest;
								L[i] = group_nearest_index;
								group[i] = gi;
								upb_at_last[i] = group_nearest;
							}
							else{
								glowers_at_last[i*ngroups + gi] = group_nearest;
							}
						}					
						else{
							arrutilv2::set_argminmin2(groupsizes[gi], distances.get(), L[i], upb_at_last[i], glowers_at_last[i*ngroups + gi]);
							L[i] += groupparts[gi];
						}
					}
					globallowers_at_last[i] = std::min(globallowers_at_last[i], glowers_at_last[i*ngroups + gi]);  // - u_delta_G[ngroups*tau[i*ngroups + gi] + gi]
				}
				//tau_globallowers[i] = modround;
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
}




template <typename TInt, typename TTau, typename TFloat>
class P21V4 : public kmeans::BaseYYMNS<TInt, TTau, TFloat>{
	
	private:
		std::unique_ptr<TTau []> tau_ub;

	protected:	
		
		TTau * const get_tau_ub(){
			return tau_ub.get();
		}		
		
		std::function<void(TInt)> update_L_glowers_upb_S_H_21v4_ati(){

			return [this](TInt ti){

				TInt x0 = (ti*this->getndata())/this->getnthreads();

				TInt modround = this->round % this->get_cutperiod();		
		
				if (modround  != 0){
					this->pll_principal_X(update_L_glowers_upb_S_H_21v4_not_modround<TInt, TTau, TFloat>, ti, this->get_u_delta_C(), this->get_u_delta_G(), this->get_u_delta_global(), this->get_C_hist(), this->get_C_l22s_hist(), this->get_ngroups(), this->get_groupparts(), this->get_groupsizes(), modround, this->get_group() + x0,this->get_L() + x0,this->get_tau() + x0*this->get_ngroups(), this->get_tau_globallowers() + x0, this->get_glowers_at_last() + x0*this->get_ngroups(), this->get_globallowers_at_last() + x0, this->get_upb_at_last() + x0, this->get_tau_ub() + x0);
				}
				
				else{
					this->pll_principal_X(update_L_glowers_upb_S_H_21v4_modround<TInt, TTau, TFloat>, ti, this->get_u_delta_C(), this->get_u_delta_G(), this->get_u_delta_global(), this->get_C_hist(), this->get_C_l22s_hist(), this->get_ngroups(), this->get_groupparts(), this->get_groupsizes(), modround, this->get_cutperiod(), this->get_group() + x0,this->get_L() + x0,this->get_tau() + x0*this->get_ngroups(), this->get_tau_globallowers() + x0,   this->get_glowers_at_last() + x0*this->get_ngroups(), this->get_globallowers_at_last() + x0, this->get_upb_at_last() + x0, this->get_tau_ub() + x0);
				}
			};				
		}
		
		/* based on arrutilv2::update_u_delta_C_from_C_C_hist_ati, just added update of u_delta_G */
		std::function<void(TInt)> update_u_delta_CGglobal_from_C_C_hist_ati(){
			return [this](TInt ti){
								
				TInt nrounds_to_proc = 1 + (this->round - 1)%this->get_cutperiod();
				TInt r0 = (ti*nrounds_to_proc)/this->getnthreads();
				TInt r1 = ((ti + 1)*nrounds_to_proc)/this->getnthreads();
				TInt local_ndcalcs = 0;
								
				arrutilv2::update_u_deltaC_from_C_C_hist(this->getncentroids(), this->getdimension(), this->get_C(), r1 - r0, this->get_C_hist() + r0*this->getncentroids()*this->getdimension(), this->get_u_delta_C() + r0*this->getncentroids(), local_ndcalcs);

				//arrutilv2::update_u_deltaC_from_C_C_hist(this->getncentroids(), this->getdimension(), this->get_C(), r1 - r0, this->get_C_hist() + r0*this->getncentroids()*this->getdimension(), this->get_C_l22s(), this->get_C_l22s_hist() + r0*this->getncentroids(), this->get_u_delta_C() + r0*this->getncentroids(), local_ndcalcs);
							
				arrutilv2::update_u_delta_G_global_from_u_delta_C(r1 - r0, this->getncentroids(), this->get_u_delta_C() + r0*this->getncentroids(), this->get_ngroups(), this->get_groupparts(), this->get_u_delta_G() + this->get_ngroups()*r0, this->get_u_delta_global() + r0);
				
				
				
				 
				if (ti == this->getnthreads() - 1){
					std::fill_n(this->get_u_delta_C() + r1*this->getncentroids(), this->getncentroids(), static_cast<TInt>(0));
					std::fill_n(this->get_u_delta_G() + r1*this->get_ngroups(), this->get_ngroups(), static_cast<TInt>(0));
					*(this->get_u_delta_global() + r1) = static_cast<TInt>(0);
				
				}
				this->ndcalcs_notX += local_ndcalcs;
			};
		}
		
		
	public:
		typedef kmeans::BaseYYMNS<TInt, TTau, TFloat> BC;
		template<typename... Args>
		P21V4(Args&&... args): BC(std::forward<Args>(args)...),
		tau_ub { new TTau [this->getndata()] }
				
		{
			this->setalgname("p21v4");
			std::fill_n(tau_ub.get(), this->getndata()*1, 0);
		}
		
		virtual TInt get_approximate_memory_requirement(){
			return BC::get_approximate_memory_requirement() + sizeof(TTau)*this->getndata();
		}
				
		virtual ~P21V4(){}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to 21v4..\n\n";
		}

		virtual void set_initialisation_tasks(){
			this->yinyang_MNS_initialisation_tasks();
		}
		
		virtual void set_C_tasks(){
			
			this->C_tasks = {

				//first task: copy C and C_l22s into hist, update C, C_l22s (split by centroids)
				arrutilv2::update_C_C_hist_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_C_hist(), this->get_C_l22s_hist(), this->round, this->get_cutperiod()),
	
				//second task: update u_delta_C and u_delta_G and u_delta_global (split by rounds)
				//(should this be a MNS base function ?)
				this->update_u_delta_CGglobal_from_C_C_hist_ati(),
				
			};
		}
		
		virtual void set_X_tasks(){
			
			this->X_tasks = {
				this->update_L_glowers_upb_S_H_21v4_ati()	
			};
		}
};

}

#endif

