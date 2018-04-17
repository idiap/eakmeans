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

#ifndef PLL_PLLYINYANGSMNBASEKMEANS_H
#define PLL_PLLYINYANGSMNBASEKMEANS_H

#include "baseYY.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{


template <typename TInt, typename TFloat>
/* "sum max norm" for lower bounds versions inherit from here */ 
class BaseYYSMN : public kmeans::BaseYY<TInt, TFloat>{
	
	private:
		std::unique_ptr<TFloat []> delta_C;
		std::unique_ptr<TFloat []> delta_G;

	protected:	
		TFloat * const get_glowers(){
			return this->get_glowers_base();
		}
		
		TFloat * const get_upb(){
			return this->get_upb_base();
		}
		
		TFloat * const get_delta_C(){
			return delta_C.get();
		}
		
		TFloat * const get_delta_G(){
			return delta_G.get();
		}
		

		
	public:
		typedef kmeans::BaseYY<TInt, TFloat> YYB;
		template<typename... Args>
		BaseYYSMN(Args&&... args): YYB(std::forward<Args>(args)...), 
		delta_C{ new TFloat [this->getncentroids()] },
		delta_G{ new TFloat [this->get_ngroups()] }
		
		{
			this->setalgname("BaseYYSMN");
		}
		
		virtual ~BaseYYSMN(){}


		virtual TInt get_approximate_memory_requirement(){
			return YYB::get_approximate_memory_requirement() + 
			sizeof(TFloat)*(
				this->getncentroids() + //delta_C
				this->get_ngroups()); //delta_G
		}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to BaseYYSMN..\n\n";
		}

		
		
		
		virtual void set_initialisation_tasks() = 0;

		
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
};

}

#endif

