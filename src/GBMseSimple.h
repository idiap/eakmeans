/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef GBMSESIMPLE_H
#define GBMSESIMPLE_H

#include "BaseGrowBatchMse.h"

namespace kmeans{
	
template <typename TInt, typename TFloat>
 
class GBMseSimple : public kmeans::BaseGrowBatchMse<TInt, TFloat>{
	
	private:		

		virtual void update_already_used(TInt x0, TInt x1, TInt ti) override final{
			arrutilv2::update_L_dn_S_H_batch(x1 - x0, this->maxpermultiplyblock, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dn() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids, this->get_sums(), this->get_counts(), this->nchanges, this->work_mutex);			
			this->ndcalcs_X += (x1 - x0)*this->ncentroids;
		}
		
		virtual void update_unused(TInt x0, TInt x1, TInt ti) override final{
			arrutilv2::update_L_dn_S_H_batch_increment_only(x1 - x0, this->maxpermultiplyblock, this->dimension, this->data + this->dimension*x0, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dn() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids, this->get_sums(), this->get_counts(), this->nchanges, this->work_mutex);
			this->ndcalcs_X += (x1 - x0)*this->ncentroids;
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
		GBMseSimple(TInt batchsize0, Args&&... args): kmeans::BaseGrowBatchMse<TInt, TFloat> (batchsize0, std::forward<Args>(args)...)		
		{
			this->algname = "GBMse Simple Dense"; 
		}
			
		virtual ~GBMseSimple(){};

};


}

#endif
