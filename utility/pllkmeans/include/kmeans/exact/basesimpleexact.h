/*
EAKMeans is a fast Exact K-means library written in C++ with 
command-line interface, shared library + header files and 
Python bindings

Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

This file is part of EAKMeans.

EAKMeans is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

EAKMeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.



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
	
			auto init_tasks_A = kmeans::BaseExact<TInt, TFloat>::makeset_C_C_l22s_L_inds0_mati();
	
			auto init_task_B = kmeans::BaseExact<TInt, TFloat>::set_S_H_ati();
			
			this->initialisation_tasks = std::move(init_tasks_A);
			this->initialisation_tasks.push_back(std::move(init_task_B));
			

			this->initialisation_tasks.push_back(
				[this](TInt ti) { std::cout << "in basesimple, this->counts[7] : " <<  this->counts[7] << std::endl; }
			);
			this->initialisation_tasks.push_back(
				[this](TInt ti) { std::cout << "in basesimple, this->sums[7] : " <<  this->sums[7] << std::endl; }
			);
			


			
		}
		
		virtual void set_X_tasks() = 0;


		virtual void set_C_tasks(){
		
			this->C_tasks = {

				arrutilv2::update_C_C_l22s_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s())

			};
		}
		
};

}

#endif

