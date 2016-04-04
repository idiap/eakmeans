#ifndef PLL_SPARSESTANDARDMINIBATCHKMEANS_H
#define PLL_SPARSESTANDARDMINIBATCHKMEANS_H

#include "basesparseminibatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class SparseStandardMiniBatch : public kmeans::BaseSparseMiniBatch<TInt, TFloat>{
	
	private:
				
		virtual void post_L_adjust_S_H() override final{
			//just update S and H by adding data which changed
			
			//TODO : these could be class variables as same calculated here as in update_L_label_changes.
			TInt data0 = this->mba.batchsize*(this->round%this->mba.nsubrounds);
			TInt data1 = std::min(data0 + this->mba.batchsize, this->ndata);
			
			sparse::increment_S_H(data0, data1, *this->ptrdata, this->get_L(), this->get_sums(), this->get_counts());
						
		}
		
	public:
					
		
		template<typename... Args>
		SparseStandardMiniBatch(Args&&... args): kmeans::BaseSparseMiniBatch<TInt, TFloat> (std::forward<Args>(args)...){
			
			this->algname = "sparse standard mini batch";
		}
		 		
		virtual ~SparseStandardMiniBatch(){};

};

}

//extern template class kmeans::SparseStandardMiniBatch<size_t, double>;
//extern template class kmeans::SparseStandardMiniBatch<size_t, float>;

#endif

