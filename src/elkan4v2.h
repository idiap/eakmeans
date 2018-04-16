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

#ifndef PLL_ELKANKMEANS_4V2_H
#define PLL_ELKANKMEANS_4V2_H

#include "baseelkan.h"


namespace kmeans{
	


template <typename TInt, typename TTau, typename TFloat>
inline void update_L_lut_S_H_4v2_not_modround(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const u_deltaC, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt * const L, TTau * const tau, TFloat * const lower_at_last, TFloat * const upper_at_last, const TInt & modround){
	
	std::unique_ptr<TFloat []> lowersx (new TFloat [ncentroids]);
	for (TInt i = 0; i < ndata; ++i){
		TFloat upperbound = upper_at_last[i] + u_deltaC[tau[i*ncentroids + L[i]]*ncentroids + L[i]];
		TInt label_before = L[i];	
		TInt ci = 0;
		while (ci < ncentroids){ //loop while upperbound is good
			lowersx[ci] = lower_at_last[i*ncentroids + ci] - u_deltaC[tau[i*ncentroids + ci]*ncentroids + ci];
			if (upperbound > lowersx[ci] && ci != L[i]){
				arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], lowersx[L[i]], ndcalcs);
				upper_at_last[i] = lowersx[L[i]];
				lower_at_last[i*ncentroids + L[i]] = lowersx[L[i]];
				tau[i*ncentroids + L[i]] = modround;
				if (upper_at_last[i] > lowersx[ci]){
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
			if (upper_at_last[i] > lowersx[ci]){
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
	


template <typename TInt, typename TTau, typename TFloat>
inline void update_L_lut_S_H_4v2_modround(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const u_deltaC, const TFloat * const C_hist, const TFloat * const C_l22s_hist, TInt * const L, TTau * const tau, TFloat * const lower_at_last, TFloat * const upper_at_last, const TInt & modround){

	for (TInt i = 0; i < ndata; ++i){
		for (TInt ci = 0; ci < ncentroids; ++ci){
			lower_at_last[i*ncentroids + ci] -= u_deltaC[tau[i*ncentroids + ci]*ncentroids + ci];
		}
		upper_at_last[i] += u_deltaC[tau[i*ncentroids + L[i]]*ncentroids + L[i]];
		std::fill_n(tau + i*ncentroids, ncentroids, static_cast<TTau>(0));
		TFloat upperbound = upper_at_last[i];
		TInt label_before = L[i];	
		TInt ci = 0;
		
		while (ci < ncentroids){
			if (upperbound > lower_at_last[i*ncentroids + ci] && ci != L[i]){
				arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], lower_at_last[i*ncentroids + L[i]], ndcalcs);
				upper_at_last[i] = lower_at_last[i*ncentroids + L[i]];
				tau[i*ncentroids + L[i]] = modround; //nec? TODO: remove (came across while constructin 21v3. ditto 6v0
				if (upper_at_last[i] > lower_at_last[i*ncentroids + ci]){
					arrutilv2::set_l2(dimension, data + i*dimension, C + ci*dimension, data_l22s[i], C_l22s[ci], lower_at_last[i*ncentroids + ci], ndcalcs);
					tau[i*ncentroids + ci] = modround; //nec? TODO: remove (came across while constructin 21v3. ditto 6v0
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
			if (upper_at_last[i] > lower_at_last[i*ncentroids + ci]){
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





/* discrepency with A4V0 with corresponing features is due to mindepend (cannot set mindepend to -1, but commenting out code I show that the ndcalcs can match exactly. */
template <typename TInt, typename TTau, typename TFloat>
class P4V2 : public kmeans::BaseElkan<TInt, TFloat>{
	
	public:
		typedef kmeans::BaseElkan<TInt, TFloat> EB;
		template<typename... Args>
		P4V2(Args&&... args): EB(std::forward<Args>(args)...), 
		
		tau { new TTau [this->getndata()*this->getncentroids()] }, 
		
		/* ignoring constants, M = Nd + Kcd + NK 
		 * where M is memory use, N is ndata, d is dimension, K is ncentroids, c is cutperiod
		 * to ensure that NS is not much heavier on memory than is SN, we must have
		 * Kcd < max (Nd, NK) --> c < max (N/K, N/d) --> c < N max (1/K, 1/d) --> c < N / min (K, d)
		 * we add the further bound that 1 <= c <= 50, resulting in the expression which follows.
		 * */
		cutperiod { std::max(static_cast<TInt>(1), ( //step 3 : make greater than or equal to 1
		std::min(static_cast<TInt>(50),  //step 2: make less than or equal to 50
		static_cast<TInt>(static_cast<TFloat>(this->getndata())/ //step 1: set to N / min (K, d)
		static_cast<TFloat>(std::min(this->getncentroids(), this->getdimension()))
		)))) },
		
		//std::min(static_cast<TInt>(50), static_cast<TInt>(1*static_cast<TFloat>(this->getndata())/static_cast<TFloat>(std::max(this->getdimension(), this->getncentroids()))))),
		
		C_hist { new TFloat [this->getncentroids()*this->getdimension()*(cutperiod)] },
		C_l22s_hist { new TFloat [this->getncentroids()*(cutperiod)] },
		u_deltaC { new TFloat [this->getncentroids()*(cutperiod + 1)] }
		
		{
			this->setalgname("p4v2");
			std::fill_n(tau.get(), this->getndata()*this->getncentroids(), 0);
			std::fill_n(u_deltaC.get(), this->getncentroids()*(cutperiod + 1), 0);
			if (this->get_cout_verbosity() > 0){
				std::cout << "cutperiod set to " << cutperiod << std::endl;
				std::cout << "memory request " << this->get_approximate_memory_requirement() / (1024.*1024. *1024.) << " GB " << std::endl;
			}


		}
		
		TInt get_cutperiod(){
			return cutperiod;
		}

		
		virtual ~P4V2(){}
		
		virtual TInt get_approximate_memory_requirement(){
			return EB::get_approximate_memory_requirement() + 
			sizeof(TTau)*this->getndata()*this->getncentroids() + //tau
			sizeof(TFloat)*cutperiod*this->getncentroids()*(this->getdimension() + 1) + // C_hist, C_hist_l22s  
			sizeof(TFloat)*this->getncentroids()*(cutperiod + 1); // u_deltaC
		}
	
	private: 
	
		std::unique_ptr<TTau []> tau;	
		TInt cutperiod;
		std::unique_ptr<TFloat []>  C_hist;
		std::unique_ptr<TFloat []>  C_l22s_hist; 
		std::unique_ptr<TFloat []>  u_deltaC;
			
	protected:
	
	
		std::function<void(TInt)> update_L_lut_S_H_4v2_ati(){

			return [this](TInt ti){

				TInt x0 = (ti*this->getndata())/this->getnthreads();

				TInt modround = this->round % this->get_cutperiod();		

				if (modround  != 0){
					this->pll_principal_X(update_L_lut_S_H_4v2_not_modround<TInt, TTau, TFloat>, ti, this->get_u_deltaC(), this->get_C_hist(), this->get_C_l22s_hist(), this->get_L() + x0,  this->get_tau() + x0*this->getncentroids(), this->get_lower_at_last() + x0*this->getncentroids(), this->get_upper_at_last() + x0, modround);
				}
				
				else{
					this->pll_principal_X(update_L_lut_S_H_4v2_modround<TInt, TTau, TFloat>, ti, this->get_u_deltaC(), this->get_C_hist(), this->get_C_l22s_hist(), this->get_L() + x0,  this->get_tau() + x0*this->getncentroids(), this->get_lower_at_last() + x0*this->getncentroids(), this->get_upper_at_last() + x0, modround);
				}
			};				
		}
	
	
		
		virtual void set_initialisation_tasks(){
			this->initialisation_tasks = this->exact_makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati();
		}
		
		virtual void set_C_tasks(){
			this->C_tasks = {
				//first task: copy C and C_l22s into hist, update C, C_l22s (split by centroids)
				arrutilv2::update_C_C_hist_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_C_hist(), this->get_C_l22s_hist(), this->round, this->get_cutperiod()),
				
				//second task: update u_deltaC (split by rounds)			
				arrutilv2::update_u_deltaC_from_C_C_hist_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_hist(), this->get_u_deltaC(), this->round, this->get_cutperiod(), this->ndcalcs_notX)
			};
		}
	
		virtual void set_X_tasks(){
			this->X_tasks = {
				//third task: update labels					
				this->update_L_lut_S_H_4v2_ati()
			};
		}
		
		
		TTau * const get_tau(){
			return tau.get();
		}
		
		
		TFloat * const get_lower_at_last(){
			//return this->get_lowers_base();
			return this->elkan_lowers_base.get();
		}
		
		TFloat * const get_upper_at_last(){
			return this->elkan_upper_base.get();
		}
		
		TFloat * get_C_hist(){
			return C_hist.get();
		}
		
		TFloat * get_C_l22s_hist(){
			return C_l22s_hist.get();
		}
		
		TFloat * const get_u_deltaC(){
			return u_deltaC.get();
		}
	

		
};

}

#endif


