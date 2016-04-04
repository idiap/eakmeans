#ifndef PLL_MINIBATCHKMEANS_H
#define PLL_MINIBATCHKMEANS_H

#include "basesimpleminibatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
//Like D Sculley, but instead of just adding newly labeled data to centroids, if the data has already been used first remove it from the centroid it was assigned to previously. This breaks the 1/t convergence to the local minimum 
class MiniBatch : public kmeans::BaseSimpleMiniBatch<TInt, TFloat>{
	
	private:
	
		virtual void update_L_S_H(TInt x0, TInt x1, TInt ti) override final{
			

			if (this->round < this->mba.nsubrounds){
				this->update_L_S_H_batch_increment_only(x0, x1, ti);
			}
			
			else{
				this->update_L_S_H_batch(x0, x1, ti);
			}
		}
		
	public:
		
		
		template<typename... Args>
		MiniBatch(Args&&... args): kmeans::BaseSimpleMiniBatch<TInt, TFloat> (std::forward<Args>(args)...) {
			this->setalgname("(Improved) Mini Batch Kmeans");
		}		
		

		virtual ~MiniBatch(){};

};


}


#endif



