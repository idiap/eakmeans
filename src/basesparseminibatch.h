/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_BASESPARSEMINIBATCHKMEANS_H
#define PLL_BASESPARSEMINIBATCHKMEANS_H

#include "basesparsekmeans.h"
#include "minibatchapp.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class BaseSparseMiniBatch : public kmeans::BaseSparseKmeans<TInt, TFloat>{
	
	private:



		//different versions for sparsestandardminibatch and sparseminibatch (my version, where not just a naive add)
		virtual void post_L_adjust_S_H() = 0;
		
		//update L, label_changes in pll on batch specified by round.
		virtual std::function<void(TInt)> update_L_label_changes_ati(){
			return [this](TInt ti){
				
				
				TInt data0 = this->mba.batchsize*(this->round%this->mba.nsubrounds);
				TInt data1 = std::min(data0 + this->mba.batchsize, this->ndata);
				TInt ndata_batch = data1 - data0;
				TInt x0 = data0 + (ti*ndata_batch)/this->nthreads;
				TInt x1 = data0 + ((ti+1)*ndata_batch)/this->nthreads;
				
				//std::cout << "\nupdating in [ " << x0 << ", " << x1 << " ] " << std::endl;
				
				
				this->where_label_changes[ti].clear(); //index, old, new.
				sparse::update_L(*this->ptrdata, x0, x1, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->where_label_changes[ti]);
				this->ndcalcs_X += this->ncentroids*(x1 - x0);
				
				std::lock_guard<std::mutex> gluk(this->work_mutex);
				this->mba.nchanges_on_batch[this->mba.subround] += this->where_label_changes[ti].size();
			};
		}
	

		virtual void set_mse() override final {
			this->minibatch_set_mse(this->mba);
		}
	
		virtual void set_summaries() override final {
			this->set_summaries_minibatch(this->mba);
		}
	
	protected:
		
		minibatchapp::MiniBatchApp<TInt> mba;		
				
		virtual void set_C_tasks() override final {
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s())
			};
		}
		
		//set S, H from first batch
		std::function<void(TInt)> set_S_H_from_initial_batch_ati(){
			return [this](TInt ti){
				if (ti == 0){
					this->set_S_H(static_cast<TInt>(0), this->mba.initialising_batch_size);
				}
			};
		}
		

		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_mati(){
			return this->minibatch_makeset_C_C_l22s_L_inds0_mati(this->mba);
		}

		virtual void set_initialisation_tasks(){
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_mati();
			auto init_task_B = this->set_S_H_from_initial_batch_ati();
			this->initialisation_tasks = std::move(init_tasks_A);
			
			this->initialisation_tasks.push_back([this](TInt ti)
			{

			});
			
			this->initialisation_tasks.push_back(std::move(init_task_B));
		}
		
		

		
		virtual void set_X_tasks(){
			this->X_tasks = {
				this->update_L_label_changes_ati(),
				[this](TInt ti){
					if (ti == 0){
						this->post_L_adjust_S_H();
					}
				},
				this->minibatch_subround_update(this->mba)
			};
		}	
	
		
	public:
		void constructor_helper(const TInt & batchsize){
	
			
			this->mba = minibatchapp::MiniBatchApp<TInt>(batchsize, this->ndata);			
			this->setalgname("Base Sparse Mini Batch Kmeans");
				
		}
		
		template<typename... Args>
		BaseSparseMiniBatch(TInt batchsize, Args&&... args): kmeans::BaseSparseKmeans<TInt, TFloat> (std::forward<Args>(args)...){
			
			this->constructor_helper(batchsize);
		}
		 		
		virtual ~BaseSparseMiniBatch(){};

};


}


#endif

