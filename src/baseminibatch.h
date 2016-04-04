#ifndef PLL_BASEMINIBATCHKMEANS_H
#define PLL_BASEMINIBATCHKMEANS_H

#include "basekmeans.h"
#include "minibatchapp.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class BaseMiniBatch : public kmeans::BaseKmeans<TInt, TFloat>{
	
	private:
	
		virtual void set_mse() override final {
			this->minibatch_set_mse(this->mba);
		}
	
		virtual void set_summaries() override final {
			this->set_summaries_minibatch(this->mba);
		}
	
	protected:
		
		minibatchapp::MiniBatchApp<TInt> mba;		
		
		TInt maxpermultiplyblock;



		virtual void set_C_tasks() = 0;

		
		//set S, H from first batch
		std::function<void(TInt)> set_S_H_ati(){
			return this->base_set_S_H_ati(static_cast<TInt> (0), this->mba.initialising_batch_size);
		}
		
		
		//Not as code reducing as the baseexact version,  but easier to understand
		template <typename Function, typename... Args>
		void mb_pll_principal_X(const Function & X_updater, TInt ti, Args&&... args){
	
			arrutilv2::pll_update_L_etc(
			//The compulsory parameters to pll_update_L_etc,
			X_updater, 
			this->ncentroids, this->dimension, this->get_sums(), this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_counts(), this->get_dcounts() + ti*this->ncentroids, this->mba.nchanges_on_batch[(this->mba.subround + 1)%this->mba.nsubrounds], this->ndcalcs_X, this->work_mutex
			//The additional parameters to pll_update_L_etc with correct offset
			, std::forward<Args>(args)...);
		}

	
		
	public:
		void constructor_helper(const TInt & batchsize){
			
			this->mba = minibatchapp::MiniBatchApp<TInt>(batchsize, this->ndata);
			this->setalgname("Base Mini Batch Kmeans");
				
			this->maxpermultiplyblock = //10000000; 
			std::max(static_cast<TInt> (1),
			static_cast<TInt> ((this->getndata() * this->getdimension())/(2 * this->getncentroids() * this->nthreads)));
	
			

	
			//TODO : move to summaries:
			std::cout << "batchsize : " << batchsize << "  nsubrounds : " << this->mba.nsubrounds << "  lastbatchsize : " << this->mba.lastbatchsize << "   maxpermultiply : " << this->maxpermultiplyblock << "  initialising_batch_size : " << this->mba.initialising_batch_size << std::endl;

		
		}
					
		/* overly hungry, consult Meyers to see how I can prevent this. 
		 * if non-standard constructor args, use variadic args ala Eli Bendersky. 
		 * Note that these won't be initialised by extern template class, so changes here will require full remake
		 * */ 
		 template<typename... Args>
		 BaseMiniBatch(TInt batchsize, Args&&... args): kmeans::BaseKmeans<TInt, TFloat> (std::forward<Args>(args)...){
			 this->constructor_helper(batchsize);
		 }
		 		
		virtual ~BaseMiniBatch(){};

};


}


#endif

