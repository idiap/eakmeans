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

#ifndef PLL_SIMPLEKMEANS_H
#define PLL_SIMPLEKMEANS_H

#include "basesimpleexact.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class SimpleKmeans1 : public kmeans::BaseSimpleExactKmeans<TInt, TFloat>{


	protected: 
	
			virtual void set_X_tasks(){
		
			this->X_tasks = {
				
				arrutilv2::update_L_S_H_ati(this->nthreads, this->ndata, this->dimension, this->data, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dsums(), this->get_dcounts(), this->get_sums(), this->get_counts(), this->nchanges, this->work_mutex, this->ndcalcs_X),
				
				
			
			};
			
		}

	public:		
		template<typename... Args>
		SimpleKmeans1(Args&&... args): kmeans::BaseSimpleExactKmeans<TInt, TFloat> (std::forward<Args>(args)...) {
			this->setalgname("simple kmeans");
		}		
		virtual ~SimpleKmeans1(){};
		
};

}

//extern template class kmeans::SimpleKmeans1<size_t, double>;
//extern template class kmeans::SimpleKmeans1<size_t, float>;

#endif
