#ifndef PLL_BASEELKANMINIBATCHKMEANS_H
#define PLL_BASEELKANMINIBATCHKMEANS_H

#include "baseminibatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class BaseElkanMiniBatch : public kmeans::BaseMiniBatch<TInt, TFloat>{

		
	protected:
			
		
	public:
		typedef kmeans::BaseMiniBatch<TInt, TFloat> BC;
		template<typename... Args>
		BaseElkanMiniBatch(Args&&... args): BC(std::forward<Args>(args)...)
		
		{
			this->assignmemory_elkan_upper_lowers();
			this->setalgname("elkan minibatch base");
		}
		
		virtual ~BaseElkanMiniBatch(){}

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



