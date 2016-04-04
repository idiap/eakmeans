/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_HAMERLYKMEANS_11V0_H
#define PLL_HAMERLYKMEANS_11V0_H

#include "basehamerly.h"

namespace kmeans{





template <typename TInt, typename TFloat>
void update_L_lower_upper_S_H_11v0(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C,  const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const delta_C,   TInt * const L, TFloat * const lower, TFloat * const upper,  const TInt & round){
	
	nchanges = 0;
	ndcalcs = 0; 


	TFloat m;
	TInt oldlabel;
	std::unique_ptr<TFloat []> distances (new TFloat [ncentroids]);
	
	//TODO: Hamerly checks that the label of data point is not max-mover (if fail, add second biggest budge). 
	TFloat max_deltaC_previous_round;
	TInt index_max_deltaC_previous_round;
	arrutilv2::set_argmaxmax(ncentroids, delta_C, index_max_deltaC_previous_round, max_deltaC_previous_round);
	for (TInt i = 0; i < ndata; ++i){
		lower[i] -= max_deltaC_previous_round;
		upper[i] += max_deltaC_previous_round;	
		m = std::max(halfminCC[L[i]], lower[i]);
		if (upper[i] > m){
			arrutilv2::set_l2(dimension, data + i*dimension,  C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upper[i], ndcalcs);		
			if (upper[i] > m){
				oldlabel = L[i];
				arrutilv2::set_rl2s(dimension, data + i*dimension, ncentroids, C, data_l22s[i], C_l22s, distances.get(), ndcalcs);
				arrutilv2::set_argminmin2nocheck(ncentroids, distances.get(), L[i], upper[i], lower[i]);
				if (L[i] != oldlabel){
					++nchanges;
					++H[L[i]];
					--H[oldlabel];
					arrutilv2::addto(dimension, data + i*dimension, S + dimension*L[i]);
					arrutilv2::subtractfrom(dimension, data + i*dimension, S + dimension*oldlabel);
				}
			}
		}
	}
}




/* discrepency in ndcalcs with a11v0 is due to a11v0 performing CC computation before first round */

template <typename TInt, typename TFloat>
class P11V0 : public kmeans::BaseHamerly<TInt, TFloat>{
	
	private:
		
			
	protected:
	
		TFloat * const get_lower(){
			return this->get_lower_base();
		}
		
		TFloat * const get_upper(){
			return this->get_upper_base();
		}
			
		std::function<void(TInt)> update_L_lower_upper_S_H_11v0_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				
				this->pll_principal_X(update_L_lower_upper_S_H_11v0<TInt, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_delta_C(), this->get_L() + x0,  this->get_lower() + x0, this->get_upper() + x0, this->round);
			};
		}
		
	public:
		typedef kmeans::BaseHamerly<TInt, TFloat> BH;
		template<typename... Args>
		P11V0(Args&&... args): BH(std::forward<Args>(args)...)
		
		{
			this->setalgname("p11v0");
		}
		virtual ~P11V0(){}

		virtual void verbose_write_additional(){
			BH::verbose_write_additional();
			/* anything else to print ? */
		}
		

		virtual void set_initialisation_tasks(){
			
			
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_lower_upper_S_H_mati();
				
			/* note : CC and halfminCC don't need to be set at this point. C and C_l22s must be set so that S&H can be set above. Maybe with smarter initialisations or initialisations with kmeans++ */
			
			}
		
		virtual void set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX),

				arrutilv2::update_CC_halfminCC_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_l22s(), this->get_CC(), this->get_halfminCC(), this->ndcalcs_notX)
			};
		}
		
		virtual void set_X_tasks(){
					
			this->X_tasks = {

			
			this->update_L_lower_upper_S_H_11v0_ati()	
			};
		}
};

}

#endif

