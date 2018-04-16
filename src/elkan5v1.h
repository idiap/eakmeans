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

#ifndef PLL_ELKANKMEANS_5V1_H
#define PLL_ELKANKMEANS_5V1_H

#include "elkan3v0.h"

namespace kmeans{


template <typename TInt, typename TFloat>
void update_L_lowers_upbs_S_H_5v1(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const delta_C, TInt * const L, TFloat * const lowers, TFloat * const upbs,  const TInt & round){
	

	nchanges = 0;
	ndcalcs = 0; 

	/* experiments with this in or out the loop show that it makes little difference (@test1) */
	arrutilv2::rank1rowupdate(ncentroids, delta_C, static_cast<TFloat>(-1.), ndata, lowers);
	
	
	for (TInt i = 0; i < ndata; ++i){
		/* (@test1) */
				
		upbs[i] += delta_C[L[i]];
		if (halfminCC[L[i]] < upbs[i]){
		TInt label_before = L[i];
			TInt ci = 0;
			while (ci < ncentroids){
				if ((L[i] != ci) && (upbs[i] > lowers[i*ncentroids + ci]) && (upbs[i] > 0.5*CC[ci*ncentroids + L[i]])){
					arrutilv2::set_l2(dimension, data + i*dimension, C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upbs[i], ndcalcs);
					lowers[i*ncentroids + L[i]] = upbs[i];
					if ((upbs[i] > lowers[i*ncentroids + ci]) && (upbs[i] > 0.5*CC[ci*ncentroids + L[i]])){
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
				if ((upbs[i] > lowers[i*ncentroids + ci]) && (upbs[i] > 0.5*CC[ci*ncentroids + L[i]])){ // (L[i] != ci) && 
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



template <typename TInt, typename TFloat>
class P5V1 : public P3V0<TInt, TFloat>{



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
		
		
		std::function<void(TInt)> update_L_lowers_upbs_S_H_5v1_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				this->pll_principal_X(update_L_lowers_upbs_S_H_5v1<TInt, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_delta_C(), this->get_L() + x0,  this->get_lowers() + x0*this->getncentroids(), this->get_upbs() + x0, this->round);
			};
		}
		
	public:
		typedef kmeans::P3V0<TInt, TFloat> PC; 
		template<typename... Args>
		P5V1(Args&&... args): PC(std::forward<Args>(args)...), 
		
		CC{ new TFloat [this->getncentroids()*this->getncentroids()]  },
		halfminCC{ new TFloat [this->getncentroids()] }
		{
			this->setalgname("p5v1");
		}
		
		virtual ~P5V1(){}

		virtual void verbose_write_additional(){
			PC::verbose_write_additional();
			/* do I want to write CC as well ? */
		}

		virtual void set_initialisation_tasks(){
			

			this->initialisation_tasks = this->exact_makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati();
						
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
				this->update_L_lowers_upbs_S_H_5v1_ati()
			};
		}
};


}

#endif
