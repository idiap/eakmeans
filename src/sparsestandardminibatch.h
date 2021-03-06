/*
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <james.newling@gmail.com>
All rights reserved.

eakmeans is a library for exact and approximate k-means written in C++ and
Python. This file is part of eakmeans. See file COPYING for more details.

This file is part of eakmeans.

eakmeans is free software: you can redistribute it and/or modify
it under the terms of the 3-Clause BSD Licence. See
https://opensource.org/licenses/BSD-3-Clause for more details.

eakmeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
COPYING for more details.
*/

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

