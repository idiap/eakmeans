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

#ifndef PLL_SIMPLEKBASEKMEANS_H
#define PLL_SIMPLEKBASEKMEANS_H

#include "baseexact.h"


namespace kmeans{

//two simple versions inherit from this class : simplebatch (distances calculated in batches of data) and simple (memory light version)
template <typename TInt, typename TFloat>
class BaseSimpleExactKmeans : public kmeans::BaseExact<TInt, TFloat>{

	public:
		
		template<typename... Args>
		/* variadic args ala Eli Bendersky */
		BaseSimpleExactKmeans(Args&&... args): kmeans::BaseExact<TInt, TFloat> (std::forward<Args>(args)...) {this->setalgname("simple base");}		
		virtual ~BaseSimpleExactKmeans(){};
			
	protected:
		virtual void set_initialisation_tasks(){
	
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_mati();
	
			auto init_task_B = this->base_set_S_H_ati(static_cast<TInt>(0), this->ndata);
					
			this->initialisation_tasks = std::move(init_tasks_A);
			this->initialisation_tasks.push_back(std::move(init_task_B));			
		}
		
		virtual void set_X_tasks() = 0;

		virtual void set_C_tasks(){		
			this->C_tasks = {
				//[](TInt ti){std::cout << "C task start " << std::endl; },
				arrutilv2::update_C_C_l22s_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s())
				//[](TInt ti){std::cout << "C task end " << std::endl; }

			};
		}		
};

}

#endif

