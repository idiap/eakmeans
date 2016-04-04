#ifndef SPARSEGBMSE3V0_H
#define SPARSEGBMSE3V0_H

#include "BaseSparseGrowBatchMse.h"

#include "alg_X_selkSN.h"

namespace kmeans{
	
template <typename TInt, typename TFloat>
 
class SparseGBMse3v1 : public kmeans::BaseSparseGrowBatchMse<TInt, TFloat>{
	
	private:		

		//updates L, dn, this->mba.nchanges_on_batch[ti], ...
		virtual void sgb_update_L_etc(TInt x0, TInt x1, TInt ti){
			
			this->where_label_changes[ti].clear(); //index, old, new.

			TInt ndcalcs_local = 0;
			kmeans::sparse_update_L_lowers_upper_where_changes_3v1(this->ncentroids, this->dimension, x0, x1, *this->ptrdata, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->gba.delta_C.get(), this->where_label_changes[ti], ndcalcs_local, this->get_L() + x0, this->get_lowers() + x0*this->ncentroids, this->get_dn() + x0);
			
			std::lock_guard<std::mutex> gluk(this->work_mutex);
			this->ndcalcs_X += ndcalcs_local;
		}

		//sets L, dn, ...		
		virtual void sgb_set_L_etc(TInt x0, TInt x1, TInt ti){

			sparse::set_L_lowers_dn(*this->ptrdata, x0, x1, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_lowers(), this->get_dn());

			std::lock_guard<std::mutex> gluk(this->work_mutex);
			this->ndcalcs_X += this->ncentroids*(x1 - x0);			
		}

	

	
	protected:
	
		TFloat * const get_lowers(){
			return this->elkan_lowers_base.get();
		}
	
	
	
		virtual void set_initialisation_tasks() override final{
			auto init_tasks_A = this->bgbmse_makeset_C_C_l22s_L_lowers_dn_inds0_mati(this->gba);		
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
			
			//this->C_tasks.push_back(
				//[this](TInt ti){
					//if (ti == 0){
						//std::cout << "\n----------------------------------------------\n";
						////TInt mincount = 100000;
						//for (TInt ci = 0; ci < this->ncentroids; ++ci){
							//std::cout << this->get_counts()[ci] << " ";
							////if (mincount > this->get_counts()[ci]){
								////mincount = this->get_counts()[ci];
							////}
						//}
						//std::cout << "\n----------------------------------------------\n";

						
					//}
				//}
			//);		
		
		
		virtual void set_L_lowers_dn(TInt x0, TInt x1) override final{
			
			////TODO : do I need to increment ndcalcs_X?
			sparse::set_L_lowers_dn(*this->ptrdata, x0, x1, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_lowers(), this->get_dn());




			//arrutilv2::set_rrl2ss_argminmins<TInt, TFloat>(x1 - x0, this->dimension, 
			//this->data + x0*this->dimension, this->ncentroids, 
			//this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), 
			//this->elkan_lowers_base.get() + x0*this->ncentroids, 
			//this->get_L() + x0, 
			//this->gbmseapp.dn.get() + x0
			
			//);
		}

	public:

		template<typename... Args>
		SparseGBMse3v1(TInt batchsize0, Args&&... args): kmeans::BaseSparseGrowBatchMse<TInt, TFloat> (batchsize0, std::forward<Args>(args)...)		
		{
			this->assignmemory_elkan_lowers(); 
			this->algname = "GBMse 3v0 Sparse (turbocharged-rho)"; 
		}
			
		virtual ~SparseGBMse3v1(){};

};


}

#endif
