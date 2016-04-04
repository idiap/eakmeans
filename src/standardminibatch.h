#ifndef PLL_STANDARDMINIBATCHKMEANS_H
#define PLL_STANDARDMINIBATCHKMEANS_H

#include "basesimpleminibatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class StandardMiniBatch : public kmeans::BaseSimpleMiniBatch<TInt, TFloat>{
	
	private:
				
		virtual void update_L_S_H(TInt x0, TInt x1, TInt ti) override final{
			this->update_L_S_H_batch_increment_only(x0, x1, ti);
		}
	

		
	public:


		template<typename... Args>
		StandardMiniBatch(Args&&... args): kmeans::BaseSimpleMiniBatch<TInt, TFloat> (std::forward<Args>(args)...) {
			this->setalgname("Standard Mini Batch Kmeans");
		}		
		 		
		virtual ~StandardMiniBatch(){};

};


}

//extern template class kmeans::StandardMiniBatch<size_t, double>;
//extern template class kmeans::StandardMiniBatch<size_t, float>;

#endif


