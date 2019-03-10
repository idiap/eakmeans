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

#ifndef PLL_EXPONIONBASEKMEANS_H
#define PLL_EXPONIONBASEKMEANS_H

#include "basehamerly.h"
#include "sortutil.h"
#include <functional>

namespace kmeans{


template <typename TInt, typename TFloat>
void update_pairs_parts_indices_halvies_halfminCC(TInt ncentroids_to_update, TInt ncentroids_total,  TInt npartitions, const TFloat * const CC, 
std::pair<TFloat, TInt> * const geometricpairs_halvies, TFloat * const partitionvalues_halvies, TInt * const geometricindices, TFloat * const halfminCC){
	sort::update_pairs_parts_indices_halvies(ncentroids_to_update, ncentroids_total, npartitions, CC, geometricpairs_halvies, partitionvalues_halvies, geometricindices, std::less<std::pair<TFloat, TInt>>());
	for (TInt ci = 0; ci < ncentroids_to_update; ++ci){
		halfminCC[ci] = partitionvalues_halvies[ci*(npartitions - 1)];
	}
}


template <typename TInt, typename TFloat>
class BaseExponion : public kmeans::BaseHamerly<TInt, TFloat>{
	
	private:
		/* for the time being we ask that one refers to r12v6 or paper for details */
		TInt npartitions; 
		std::unique_ptr<std::pair<TFloat, TInt> []> geometricpairs_halvies;
		std::unique_ptr< TFloat []> partitionvalues_halvies;
		std::unique_ptr< TInt []> geometricindices;

			
	protected:
					
		std::pair<TFloat, TInt> * const get_geometricpairs_halvies(){
			return geometricpairs_halvies.get();
		}
		
		TFloat * const get_partitionvalues_halvies(){
			return partitionvalues_halvies.get();
		}
		
		TInt * const get_geometricindices(){
			return geometricindices.get();
		}
		
		
		std::function<void(TInt)> update_pairs_parts_indices_halvies_halfminCC_ati(){
			return [this](TInt ti){
				TInt c0 = (ti*this->getncentroids())/this->getnthreads();
				TInt c1 = ((ti+1)*this->getncentroids())/this->getnthreads();
				update_pairs_parts_indices_halvies_halfminCC(c1 - c0, this->getncentroids(), this->get_npartitions(), this->get_CC() + c0*this->getncentroids(), this->get_geometricpairs_halvies() + (this->getncentroids() - 1)*c0, this->get_partitionvalues_halvies() + (this->get_npartitions() - 1)*c0, this->get_geometricindices() + (this->getncentroids() - 1)*c0, this->get_halfminCC() + c0);
			};
		} 

		
		
	public:
		typedef kmeans::BaseHamerly<TInt, TFloat> BH;
		template<typename... Args>
		BaseExponion(Args&&... args): BH(std::forward<Args>(args)...),\
		
		
		npartitions { static_cast<TInt>(std::floor(std::log2(static_cast<TFloat>(this->getncentroids() - 1)))) },
		geometricpairs_halvies { new std::pair<TFloat, TInt> [this->getncentroids()*(this->getncentroids() - 1)] },
		partitionvalues_halvies { new TFloat [this->getncentroids()*(this->get_npartitions() - 1)] },
		geometricindices { new TInt [this->getncentroids()*(this->getncentroids() - 1)] }
		
		{
			this->setalgname("exponion");
			//as per algorithmhamerlyoriginals.h,  geometricpairs_halvies: each row contains (unitialised, i) for all i != row. 
			//TODO: move this to initialisation tasks in parallelised form 
			for (TInt ci = 0; ci < this->getncentroids(); ++ci){
				for (TInt cj = 0; cj < this->getncentroids(); ++ cj){
					if (cj < ci){
						geometricpairs_halvies[ci*(this->getncentroids() - 1) + cj].second = cj;
					}
					else if (cj > ci){
						geometricpairs_halvies[ci*(this->getncentroids() - 1) + cj -1].second = cj;
					}
				}
			}
		}
		virtual ~BaseExponion(){}

		virtual void verbose_write_additional(){
			BH::verbose_write_additional();
			/* anything else to print ? */
		}
		
		TInt get_npartitions(){
			return npartitions;
		}


		void exponion_set_initialisation_tasks(){
				
			
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_lower_upper_S_H_mati();
			/* note : as usual CC and halfminCC don't need to be set at this point. SImilarly, 
			 * geometricpairs_halvies has had the necessary skeleton put in place in constructor
			 * partitionvalues_halvies and geometricindices will be populated from scratch in run of C_tasks,
			 * no need to initialise here
			 * */			
		}
		
		virtual void set_initialisation_tasks() = 0;
		
		void exponion_set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX),

				arrutilv2::update_CC_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_l22s(), this->get_CC(), this->ndcalcs_notX),
				
				/* update :
				 * geometricpairs_halvies
				 * partitionvalues_halvies
				 * geometricindices
				 * and halfminCC
				 */
				 
				this->update_pairs_parts_indices_halvies_halfminCC_ati()
			};
		}
		
		virtual void set_C_tasks() = 0;
		
		virtual void set_X_tasks() = 0;
};

}

#endif

