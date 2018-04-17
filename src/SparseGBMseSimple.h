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

//blindly based on GBMseSimple, again very similar.

#ifndef SPARSEGBMSESIMPLE_H
#define SPARSEGBMSESIMPLE_H

#include "BaseSparseGrowBatchMse.h"

namespace kmeans{
	
template <typename TInt, typename TFloat>
 
class SparseGBMseSimple : public kmeans::BaseSparseGrowBatchMse<TInt, TFloat>{
	
	private:		

		//updates L, dn, this->mba.nchanges_on_batch[ti]		
		virtual void sgb_update_L_etc(TInt x0, TInt x1, TInt ti){
			
			this->where_label_changes[ti].clear(); //index, old, new.
			sparse::update_L_dn(*this->ptrdata, x0, x1, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dn(), this->where_label_changes[ti]);
			
			std::lock_guard<std::mutex> gluk(this->work_mutex);
			//this->nchanges += this->where_label_changes[ti].size();
			this->ndcalcs_X += this->ncentroids*(x1 - x0);			
		}

		//sets L, dn		
		virtual void sgb_set_L_etc(TInt x0, TInt x1, TInt ti){

			sparse::set_L_dn(*this->ptrdata, x0, x1, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dn());

			std::lock_guard<std::mutex> gluk(this->work_mutex);
				this->ndcalcs_X += this->ncentroids*(x1 - x0);			
		}

	

	
	protected:
	
		virtual void set_initialisation_tasks() override final{
			auto init_tasks_A = this->bgbmse_makeset_C_C_l22s_L_dn_inds0_mati(this->gba);
			auto init_task_B = this->base_set_S_H_ati(static_cast<TInt>(0), this->gba.ndata_active);	
			this->initialisation_tasks.insert(this->initialisation_tasks.end(), init_tasks_A.begin(), init_tasks_A.end());
			this->initialisation_tasks.push_back(init_task_B);
		}
		
		virtual void set_C_tasks() override final{
			this->C_tasks = {};
			
			
			//we use this codefrag to confirm that S and H are correctly set.
			if (true == false){
				this->C_tasks.push_back(
					[this](TInt ti){
						if (ti == 0){
							sparse::todense::set_S_H(*this->ptrdata, static_cast<TInt> (0), this->gba.ndata_active, this->ncentroids, this->get_L(), this->get_sums(), this->get_counts());
						}
					}
				);
			}



			this->C_tasks.push_back(			
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->gba.delta_C.get(), this->ndcalcs_notX)
			);
			this->C_tasks.push_back(
				this->set_mse_sse_by_cluster_ati(this->gba, this->gbmseapp)
			);				
			this->C_tasks.push_back(
				this->update_ndata_active_ati(this->gba)
			);
			
		}

	public:

		template<typename... Args>
		SparseGBMseSimple(TInt batchsize0, Args&&... args): kmeans::BaseSparseGrowBatchMse<TInt, TFloat> (batchsize0, std::forward<Args>(args)...)		
		{
			this->algname = "GBMse Simple Sparse"; 
		}
			
		virtual ~SparseGBMseSimple(){};

};


}

#endif
