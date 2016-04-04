#ifndef PLL_BASESIMPLEMINIBATCHKMEANS_H
#define PLL_BASESIMPLEMINIBATCHKMEANS_H

#include "baseminibatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class BaseSimpleMiniBatch : public kmeans::BaseMiniBatch<TInt, TFloat>{


	private:
	
		virtual void update_L_S_H(TInt x0, TInt x1, TInt ti) = 0;
		
		virtual std::function<void(TInt)> update_L_S_H_ati() override final{
			
			return [this](TInt ti){
				//the batch to use this round (same for all threads)				
				TInt ind0 = this->mba.batchsize*((this->mba.subround + 1) % this->mba.nsubrounds); //(this->round%this->mba.nsubrounds); //not this->mba.subround!
				TInt ind1 = std::min(ind0 + this->mba.batchsize, this->ndata);
				TInt thisbatchsize = ind1 - ind0;
				
				//absolute indices of data to process on this threads
				TInt x0 = ind0 + (ti*thisbatchsize)/this->nthreads;
				TInt x1 = ind0 + ((ti + 1)*thisbatchsize)/this->nthreads;
			
				this->update_L_S_H(x0, x1, ti);
				
				this->ndcalcs_X += (x1 - x0)*this->ncentroids;
			};
		}






	protected:
	
	
		std::unique_ptr<TFloat []> delta_C;	
				
		virtual void set_C_tasks() override final {
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->delta_C.get(), this->ndcalcs_notX)
			};
		}
		

		//some initialisation, using the first batch if necessary
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_mati(){
			return this->minibatch_makeset_C_C_l22s_L_inds0_mati(this->mba);
		}

		
		
		
		virtual void set_initialisation_tasks() override final{
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_mati();
			auto init_task_B = this->set_S_H_ati();
			this->initialisation_tasks = std::move(init_tasks_A);
			this->initialisation_tasks.push_back(std::move(init_task_B));
		}
		
		

		
		virtual void set_X_tasks() override final {
			this->X_tasks = {
				
				//[this](TInt ti){
					//std::cout << "\nnchanges_on_batch" << std::endl;
					//for (TInt a = 0; a < this->mba.nsubrounds; ++a){
						//std::cout << this->mba.nchanges_on_batch[a] << " ";
					//}
					//std::cout << std::endl;
				//},
				
				this->update_L_S_H_ati(),
				this->minibatch_subround_update(this->mba)
			};
		}
		
		
		void update_L_S_H_batch_increment_only(TInt x0, TInt  x1, TInt ti){
			arrutilv2::update_L_S_H_batch_increment_only(x1-x0, this->maxpermultiplyblock, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_sums(), this->get_counts(), this->mba.nchanges_on_batch[(this->mba.subround + 1)%this->mba.nsubrounds], this->work_mutex);
		}
		
		
		void update_L_S_H_batch(TInt x0, TInt  x1,  TInt ti){
			arrutilv2::update_L_S_H_batch(x1-x0, this->maxpermultiplyblock, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_sums(), this->get_counts(), this->mba.nchanges_on_batch[(this->mba.subround + 1)%this->mba.nsubrounds], this->work_mutex);
		}
	


		
		
	public:

	template<typename... Args>
	BaseSimpleMiniBatch(Args&&... args): kmeans::BaseMiniBatch<TInt, TFloat> (std::forward<Args>(args)...){
	 			this->delta_C = std::unique_ptr<TFloat []> (new TFloat [this->ncentroids]);
	}
		
	virtual ~BaseSimpleMiniBatch(){};
	
};
}

#endif
