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


