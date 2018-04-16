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



