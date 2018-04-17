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



