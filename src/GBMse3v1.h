#ifndef GBMse3v1_H
#define GBMse3v1_H



#include "BaseGrowBatchMse.h"
#include "alg_X_selkSN.h"
#include "arrutilv2l3.h"

namespace kmeans{
	
//The updating with 3v1 is similar to that of 3v0, but takes advantage of the fact that upper bounds are always tight. (dn is always the distance to the nearest centroid).
	
template <typename TInt, typename TFloat>
 
class GBMse3v1 : public kmeans::BaseGrowBatchMse<TInt, TFloat>{
	
	private:		

		virtual void update_already_used(TInt x0, TInt x1, TInt ti) override final{
			
			this->gb_pll_principal_X(
			kmeans::update_L_lowers_upper_S_H_3v1<TInt, TFloat>, 
			ti,
			x1 - x0,
			this->data + x0*this->dimension,
			this->get_C(),
			this->get_data_l22s() + x0, 
			this->get_C_l22s(), 
			this->gba.delta_C.get(), 
			this->get_L() + x0, 
			this->get_lowers() + x0*this->ncentroids, 
			this->get_dn() + x0, 
			this->round);
			
		}
		
		virtual void update_unused(TInt x0, TInt x1, TInt ti) override final{
					
			arrutilv2::set_L_lowers_dn_and_increment_S_H(x1 - x0, this->dimension, this->data + this->dimension*x0, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_lowers() + x0*this->ncentroids, this->get_dn() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids, this->get_sums(), this->get_counts(), this->nchanges, this->work_mutex);
			this->ndcalcs_X += (x1 - x0)*this->ncentroids;
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
		
		virtual void set_L_lowers_dn(TInt x0, TInt x1) override final{
						
			arrutilv2::set_rrl2ss_argminmins<TInt, TFloat>(x1 - x0, this->dimension, 
			this->data + x0*this->dimension, this->ncentroids, 
			this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), 
			this->elkan_lowers_base.get() + x0*this->ncentroids, 
			this->get_L() + x0, 
			this->gbmseapp.dn.get() + x0
			
			);
		}

	public:

		template<typename... Args>
		GBMse3v1(TInt batchsize0, Args&&... args): kmeans::BaseGrowBatchMse<TInt, TFloat> (batchsize0, std::forward<Args>(args)...)		
		{
			this->assignmemory_elkan_lowers(); 
			this->algname = "GBMse Elkan 3v1";
		}
			
		virtual ~GBMse3v1(){};

};

}

#endif
