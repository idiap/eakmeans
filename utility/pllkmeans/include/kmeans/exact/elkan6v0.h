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

#ifndef PLL_ELKANKMEANS_6V0_H
#define PLL_ELKANKMEANS_6V0_H


#include "elkan4v2.h"

namespace kmeans{



template <typename TInt, typename TTau, typename TFloat>
inline void update_L_lut_S_H_6v0_not_modround(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const u_deltaC, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt * const L, TTau * const tau, TFloat * const lower_at_last, TFloat * const upper_at_last, const TInt & modround){
	
	std::unique_ptr<TFloat []> lowersx (new TFloat [ncentroids]);
	for (TInt i = 0; i < ndata; ++i){
		TFloat upperbound = upper_at_last[i] + u_deltaC[tau[i*ncentroids + L[i]]*ncentroids + L[i]];
		if (halfminCC[L[i]] < upperbound){
			TInt label_before = L[i];	
			TInt ci = 0;
			while (ci < ncentroids){ //loop while upperbound is good
				lowersx[ci] = lower_at_last[i*ncentroids + ci] - u_deltaC[tau[i*ncentroids + ci]*ncentroids + ci];
				if (upperbound > lowersx[ci] && ci != L[i] && (upperbound > 0.5*CC[ci*ncentroids + L[i]])){
					arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], lowersx[L[i]], ndcalcs);
					upper_at_last[i] = lowersx[L[i]];
					lower_at_last[i*ncentroids + L[i]] = lowersx[L[i]];
					tau[i*ncentroids + L[i]] = modround;
					if (upper_at_last[i] > lowersx[ci] && (upper_at_last[i] > 0.5*CC[ci*ncentroids + L[i]])){
						arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lowersx[ci], ndcalcs);
						lower_at_last[i*ncentroids + ci] = lowersx[ci];
						tau[i*ncentroids + ci] = modround;
						if (upper_at_last[i] > lowersx[ci]){
							upper_at_last[i] = lowersx[ci];
							L[i] = ci;
						}
					}
					++ci;
					break;
				}
				++ci;
			}
		
			while (ci < ncentroids){ //upperbound is now tight
				lowersx[ci] = lower_at_last[i*ncentroids + ci] - u_deltaC[tau[i*ncentroids + ci]*ncentroids + ci];
				if (upper_at_last[i] > lowersx[ci] && (upper_at_last[i] > 0.5*CC[ci*ncentroids + L[i]])){
					arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lowersx[ci], ndcalcs);
					lower_at_last[i*ncentroids + ci] = lowersx[ci];
					tau[i*ncentroids + ci] = modround;
					if (upper_at_last[i] > lowersx[ci]){
						upper_at_last[i] = lowersx[ci];
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
	


template <typename TInt, typename TTau, typename TFloat>	
inline void update_L_lut_S_H_6v0_modround(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C,  const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const u_deltaC, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt * const L, TTau * const tau, TFloat * const lower_at_last, TFloat * const upper_at_last, const TInt & modround){

	for (TInt i = 0; i < ndata; ++i){
		for (TInt ci = 0; ci < ncentroids; ++ci){
			lower_at_last[i*ncentroids + ci] -= u_deltaC[tau[i*ncentroids + ci]*ncentroids + ci];
		}
		upper_at_last[i] += u_deltaC[tau[i*ncentroids + L[i]]*ncentroids + L[i]];
		std::fill_n(tau + i*ncentroids, ncentroids, static_cast<TTau>(0));
		TFloat upperbound = upper_at_last[i];
		
		if (halfminCC[L[i]] < upperbound){
			TInt label_before = L[i];	
			TInt ci = 0;
			while (ci < ncentroids){
				if (upperbound > lower_at_last[i*ncentroids + ci] && ci != L[i] && (upperbound > 0.5*CC[ci*ncentroids + L[i]])){
					arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], lower_at_last[i*ncentroids + L[i]], ndcalcs);
					upper_at_last[i] = lower_at_last[i*ncentroids + L[i]];
					tau[i*ncentroids + L[i]] = modround;
					if (upper_at_last[i] > lower_at_last[i*ncentroids + ci] && (upper_at_last[i] > 0.5*CC[ci*ncentroids + L[i]])){
						arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lower_at_last[i*ncentroids + ci], ndcalcs);
						tau[i*ncentroids + ci] = modround;
						if (upper_at_last[i] > lower_at_last[i*ncentroids + ci]){
							upper_at_last[i] = lower_at_last[i*ncentroids + ci];
							L[i] = ci;
						}
					}
					++ci;
					break;
				}
				++ci;
			}
		
			while (ci < ncentroids){
				if (upper_at_last[i] > lower_at_last[i*ncentroids + ci]  && (upper_at_last[i] > 0.5*CC[ci*ncentroids + L[i]])){
					arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lower_at_last[i*ncentroids + ci], ndcalcs);
					tau[i*ncentroids + ci] = modround;
					if (upper_at_last[i] > lower_at_last[i*ncentroids + ci]){
						upper_at_last[i] = lower_at_last[i*ncentroids + ci];
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


template <typename TInt, typename TTau, typename TFloat>
void update_L_lut_S_H_6v0(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const CC,  const TFloat * const halfminCC, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const u_deltaC, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt * const L, TTau * const tau, TFloat * const lower_at_last, TFloat * const upper_at_last, const TInt & round, TInt cutperiod){
	
	TInt modround = round%cutperiod;

	if (modround != 0){
		update_L_lut_S_H_6v0_not_modround<TInt, TTau, TFloat>(ncentroids, dimension, S, H, nchanges, ndcalcs, ndata, data, C, CC, halfminCC, data_l22s, C_l22s, u_deltaC, C_hist, C_l22s_hist, L, tau, lower_at_last, upper_at_last, modround);
	}
	
	else{
		update_L_lut_S_H_6v0_modround< TInt, TTau, TFloat>(ncentroids, dimension, S, H, nchanges, ndcalcs, ndata, data, C, CC, halfminCC, data_l22s, C_l22s, u_deltaC, C_hist, C_l22s_hist, L, tau, lower_at_last, upper_at_last, modround);
	}
}






template <typename TInt, typename TTau, typename TFloat>
class P6V0 : public P4V2<TInt, TTau, TFloat>{

	private:
		std::unique_ptr<TFloat []> CC;
		std::unique_ptr<TFloat []> halfminCC;	
		
	protected:
		TFloat * const get_CC(){
			return CC.get();
		}
		
		TFloat * const get_halfminCC(){
			return halfminCC.get();
		}
		
		std::function<void(TInt)> update_L_lut_S_H_6v0_ati(){

			return [this](TInt ti){

				TInt x0 = (ti*this->getndata())/this->getnthreads();

				TInt modround = this->round % this->get_cutperiod();		

				if (modround  != 0){
					this->pll_principal_X(update_L_lut_S_H_6v0_not_modround<TInt, TTau, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_u_deltaC(), this->get_C_hist(), this->get_C_l22s_hist(), this->get_L() + x0,  this->get_tau() + x0*this->getncentroids(), this->get_lower_at_last() + x0*this->getncentroids(), this->get_upper_at_last() + x0, modround);
				}
				
				else{
					this->pll_principal_X(update_L_lut_S_H_6v0_modround<TInt, TTau, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_u_deltaC(), this->get_C_hist(), this->get_C_l22s_hist(), this->get_L() + x0,  this->get_tau() + x0*this->getncentroids(), this->get_lower_at_last() + x0*this->getncentroids(), this->get_upper_at_last() + x0, modround);
				}
			};				
		}
		
		
	public:
		typedef kmeans::P4V2<TInt, TTau, TFloat> PC; //parent class
		template<typename... Args>
		P6V0(Args&&... args): PC(std::forward<Args>(args)...), 
		
		CC{ new TFloat [this->getncentroids()*this->getncentroids()]  },
		halfminCC{ new TFloat [this->getncentroids()] }
		{
			this->setalgname("p6v0");
		}
		
		virtual ~P6V0(){}

		virtual void verbose_write_additional(){
			PC::verbose_write_additional();
			/* do I want to write CC as well ? */
		}

		virtual void set_initialisation_tasks(){
			

			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati();
				
			/* note again : CC and halfminCC don't need to be set at this point. C and C_l22s must be set so that S&H can be set above. Maybe with smarter initialisations or initialisations with kmeans++ */
			
		}
		
		virtual void set_C_tasks(){
			this->C_tasks = {		
				//first task: copy C and C_l22s into hist, update C, C_l22s (split by centroids)
				arrutilv2::update_C_C_hist_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_C_hist(), this->get_C_l22s_hist(), this->round, this->get_cutperiod()),
				
	
				//second task: update u_deltaC (split by rounds)			
				arrutilv2::update_u_deltaC_from_C_C_hist_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_hist(), this->get_u_deltaC(), this->round, this->get_cutperiod(), this->ndcalcs_notX),
				
				//third task : TODO: merge with second task.
				arrutilv2::update_CC_halfminCC_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_l22s(), this->get_CC(), this->get_halfminCC(), this->ndcalcs_notX)
			};
		}
		
		virtual void set_X_tasks(){
			this->X_tasks = {
				this->update_L_lut_S_H_6v0_ati()
			};
		}
};


}

#endif
