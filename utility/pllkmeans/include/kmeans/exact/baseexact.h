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

#ifndef PLL_BASEEXACTKMEANSTRUE_H
#define PLL_BASEEXACTKMEANSTRUE_H

#include "basekmeans.h"


namespace kmeans{

template <typename TInt, typename TFloat>
class BaseExact : public kmeans::BaseKmeans<TInt, TFloat> {
	
	private:
		virtual void set_summaries() {
			this->set_summaries_exact();
		}
		
		
		virtual void set_mse() override {
			this->mse = arrutilv2::getmeanl22at(this->ncentroids, this->dimension, this->get_C(), this->ndata, this->data, this->get_L(), this->get_C_l22s(), this->get_data_l22s());
		}

	protected:
		virtual void set_initialisation_tasks() = 0;
		virtual void set_X_tasks() = 0;
		virtual void set_C_tasks() = 0;
	
		template <typename Function, typename... Args>
		void pll_principal_X(const Function & X_updater, TInt ti, Args&&... args){
			this->base_pll_principal_X(static_cast<TInt> (0), this->ndata, X_updater, ti, std::forward<Args>(args)...);
		}
		
		std::function<void(TInt)> set_S_H_ati(){
			return this->base_set_S_H_ati(static_cast<TInt> (0), this->ndata);
		}
		
		/* return a vector of functions which take thread number as argument and do work, in this case setting C and C_l22s and this->inds0. Appropriate for paralellised framework */
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_inds0_mati(){
			return this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), this->ndata);
		}

		
		
	
	public:
		template<typename... Args>
		BaseExact(Args&&... args): kmeans::BaseKmeans<TInt, TFloat>(std::forward<Args>(args)...){
			
		}
		
		virtual ~BaseExact(){}

		
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_mati(){
		
		//will not be used by all who inherit	 (only simplest and basesimple)
		std::vector<std::function<void(TInt)> > tasks;
			tasks = this->makeset_C_C_l22s_inds0_mati();			
			tasks.push_back(
			[this](TInt ti){
				TInt local_ndcalcs = 0;
				TInt x0 = (ti*this->ndata)/this->nthreads;
				TInt x1 = ((ti+1)*this->ndata)/this->nthreads;
				arrutilv2::set_rargmins(x1 - x0, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, local_ndcalcs);
				this->ndcalcs_notX += local_ndcalcs;
				}
			);
			
			return tasks;
		}

	
};

}


//extern template class kmeans::BaseExact<size_t, double>;
//extern template class kmeans::BaseExact<size_t, float>;



#endif



