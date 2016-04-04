/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef BASEGROWBATCHMSE_H
#define BASEGROWBATCHMSE_H

#include "BaseGrowBatch.h"
#include "growbatchapp.h"

#include <algorithm>
#include <functional>
namespace kmeans{
	
template <typename TInt, typename TFloat>
/* A type of GrowBatch, so batch size doubles when determined to be appropriate.
 * Specifically, we monitor 
 * (1) mse per cluster, 
 * (2) delta_C per cluster 
 * and if *median* of mse/delta_C > *1*, double (while can) see basecluster function for details (in basecluster so that sparse can use as well)
 * */
 
class BaseGrowBatchMse : public kmeans::BaseGrowBatch<TInt, TFloat>{
	
	private:		
		virtual void set_mse() override final {
			this->gbmse_set_mse(this->gba, this->gbmseapp);
		}
		
	protected:

		growbatchapp::GBMseApp<TInt, TFloat> gbmseapp;
		
		TFloat * const get_dn(){
			return this->gbmseapp.dn.get();
		}
	
		virtual bool should_double() override final{
			return this->gbmse_should_double(this->gba, this->gbmseapp.mse_by_cluster.data());
		}
		
		virtual void set_L_dn(TInt x0, TInt x1) override final {
			TInt local_ndcalcs = 0;			
			arrutilv2::set_rargminmins(x1 - x0, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dn() + x0, local_ndcalcs);
			this->ndcalcs_X += local_ndcalcs;
		}
			
		virtual void set_initialisation_tasks() = 0;
		
		virtual void set_X_tasks() override final{
			/* using function as defined in basedensecentroidkmeans */
			this->X_tasks = this->bgbmse_update_L_dn_etc_S_H_batch_switch_mati(this->gba); 		
		}
		
		virtual void set_C_tasks() = 0;
		
		public:

			template<typename... Args>
			BaseGrowBatchMse(TInt batchsize0, Args&&... args): kmeans::BaseGrowBatch<TInt, TFloat> (batchsize0, std::forward<Args>(args)...)		
			{
				this->BGBM_constructor_helper(this->gbmseapp);			
			}
				
			virtual ~BaseGrowBatchMse(){};

};
}

#endif
