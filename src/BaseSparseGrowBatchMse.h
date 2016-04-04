//based blindly on BaseSparseGrowBatchMse.h, almost identical. 

#ifndef BASESPARSEGROWBATCHMSE_H
#define BASESPARSEGROWBATCHMSE_H

#include "BaseSparseGrowBatch.h"
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
 
class BaseSparseGrowBatchMse : public kmeans::BaseSparseGrowBatch<TInt, TFloat>{
	
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
			sparse::set_L_dn(*this->ptrdata, x0, x1, this->ncentroids, this->get_C(), 
			this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dn());
			this->ndcalcs_X += (x1 - x0)*this->ncentroids;
		}
			
		virtual void set_initialisation_tasks() = 0;
				
		virtual void set_C_tasks() = 0;
		
		public:

			template<typename... Args>
			BaseSparseGrowBatchMse(TInt batchsize0, Args&&... args): kmeans::BaseSparseGrowBatch<TInt, TFloat> (batchsize0, std::forward<Args>(args)...)		
			{
				this->BGBM_constructor_helper(this->gbmseapp);			
			}
				
			virtual ~BaseSparseGrowBatchMse(){};

};
}

#endif
