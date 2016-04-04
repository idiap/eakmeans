#ifndef PLL_BASESPARSEEXACT_H
#define PLL_BASESPARSEEXACT_H

#include "basesparsekmeans.h"
#include <stdexcept>

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseSparseExact : public kmeans::BaseSparseKmeans<TInt, TFloat> {
	
	private: 
		virtual void set_mse() override{
			TFloat sse = 0;
			for (TInt i = 0; i < this->ndata; ++i){
				sse += 
				this->data_l22s[i] + this->C_l22s[this->L[i]]				
				-2.*sparse::get_inner(this->ptrdata->starts[i+1] - this->ptrdata->starts[i], 
				this->ptrdata->indices.data() + this->ptrdata->starts[i], 
				this->ptrdata->values.data() + this->ptrdata->starts[i],
				this->get_C() + this->dimension*this->L[i]);
			}
			this->mse =  sse / static_cast<TFloat> (this->ndata);
		}

		
	public:
		template<typename... Args>
		BaseSparseExact(Args&&... args): kmeans::BaseSparseKmeans<TInt, TFloat> (std::forward<Args>(args)...) {
			this->setalgname("base-sparse-exact-kmeans");
		}
		
		virtual ~BaseSparseExact(){};
	
	protected:
		virtual void set_initialisation_tasks() = 0;
		virtual void set_X_tasks() = 0;		
		virtual void set_C_tasks() = 0;
		
		virtual void set_summaries(){
			this->set_summaries_exact();
		}
		



		//A hack as no pllsation as suggested by ati suffix
		std::function<void(TInt)> set_S_H_ati(){
			return this->base_set_S_H_ati(static_cast<TInt>(0), this->ndata);
		}

		
		virtual void verbose_write_additional(){
			throw std::runtime_error("verbose_write_additional needs implementing in basesparseexact");
		}

	

};
} 


#endif
