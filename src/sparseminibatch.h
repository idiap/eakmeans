#ifndef PLL_SPARSEMINIBATCHKMEANS_H
#define PLL_SPARSEMINIBATCHKMEANS_H

#include "basesparseminibatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class SparseMiniBatch : public kmeans::BaseSparseMiniBatch<TInt, TFloat>{
	
	private:
				
		virtual void post_L_adjust_S_H() override final{

			if (this->round < this->mba.nsubrounds){
				TInt data0 = this->mba.batchsize*(this->round%this->mba.nsubrounds);
				TInt data1 = std::min(data0 + this->mba.batchsize, this->ndata);
				sparse::increment_S_H(data0, data1, *this->ptrdata, this->get_L(), this->get_sums(), this->get_counts());
			}
			
			else{
				sparse::update_S_H_from_label_changes(*this->ptrdata, this->where_label_changes, this->get_sums(), this->get_counts());
			}						
		}
				
		
	public:
					
		
		template<typename... Args>
		SparseMiniBatch(Args&&... args): kmeans::BaseSparseMiniBatch<TInt, TFloat> (std::forward<Args>(args)...){
				this->algname = "sparse mini batch";
		}
		 		
		virtual ~SparseMiniBatch(){};

};

}

#endif
