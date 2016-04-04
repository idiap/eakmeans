#ifndef PLL_BASESPARSEELKANKMEANS__H
#define PLL_BASESPARSEELKANKMEANS__H

#include "basesparseexact.h"
#include "sparseutil.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseSparseElkan : public kmeans::BaseSparseExact<TInt, TFloat>{
		
	protected:
			
		virtual void set_upper_lowers_L(TInt x0, TInt x1) override final{ /* from basedensecentroidkmeans */
			for (TInt i = x0; i < x1; ++i){								
				sparse::set_argminmin_rl2s(this->ptrdata->starts[i+1] - this->ptrdata->starts[i], this->ptrdata->indices.data() + this->ptrdata->starts[i],  this->ptrdata->values.data() + this->ptrdata->starts[i], this->ptrdata->dimension, this->ncentroids, this->get_C(), this->data_l22s[i], this->get_C_l22s(), this->L[i], this->elkan_upper_base[i], this->elkan_lowers_base.get() + i*this->ncentroids);
			}
			
			this->ndcalcs_X += this->ncentroids*(x1 - x0);
		}
		
	public:
		typedef kmeans::BaseSparseExact<TInt, TFloat> BC;
		template<typename... Args>
		BaseSparseElkan(Args&&... args): BC(std::forward<Args>(args)...)
		
		{
			this->assignmemory_elkan_upper_lowers();
			this->setalgname("sparse elkan base");
		}
		
		virtual ~BaseSparseElkan(){}

		virtual void verbose_write_additional() override {}
		virtual void set_initialisation_tasks() = 0;
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
		
		virtual TInt get_approximate_memory_requirement(){
			return BC::get_approximate_memory_requirement() + this->get_elkan_base_memory();
		}
};

}

#endif



