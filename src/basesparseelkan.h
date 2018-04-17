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

#ifndef PLL_BASESPARSEELKANKMEANS__H
#define PLL_BASESPARSEELKANKMEANS__H

#include "basesparseexact.h"
#include "sparseutil.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseSparseElkan : public kmeans::BaseSparseExact<TInt, TFloat>{
		
	protected:
			
		virtual void set_upper_lowers_L(TInt x0, TInt x1) override final{ /* from basedensecentroidkmeans */
			for (TInt i = x0; i < x1; ++i){								
				sparse::set_argminmin_rl2s(this->ptrdata->starts[i+1] - this->ptrdata->starts[i], this->ptrdata->indices.data() + this->ptrdata->starts[i],  this->ptrdata->values.data() + this->ptrdata->starts[i], this->ptrdata->dimension, this->ncentroids, this->get_C(), this->data_l22s[i], this->get_C_l22s(), this->L[i], this->elkan_upper_base[i], this->elkan_lowers_base.get() + i*this->ncentroids);
			}
			
			this->ndcalcs_X += this->ncentroids*(x1 - x0);
		}
		
	public:
		typedef kmeans::BaseSparseExact<TInt, TFloat> BC;
		template<typename... Args>
		BaseSparseElkan(Args&&... args): BC(std::forward<Args>(args)...)
		
		{
			this->assignmemory_elkan_upper_lowers();
			this->setalgname("sparse elkan base");
		}
		
		virtual ~BaseSparseElkan(){}

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



