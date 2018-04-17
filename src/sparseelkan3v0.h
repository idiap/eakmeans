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

#ifndef PLL_SPARSEELKANKMEANS_3V0_H
#define PLL_SPARSEELKANKMEANS_3V0_H

#include "basesparseelkan.h"
#include "alg_X_selkSN.h"

namespace kmeans{

/* discrepency in ndcalcs as compared to a3v0 due to not computing CC initially (I propose) */

template <typename TInt, typename TFloat>
class SP3V0 : public kmeans::BaseSparseElkan<TInt, TFloat>{
				
	protected:
		TFloat * const get_lowers(){
			return this->elkan_lowers_base.get();
		}
		
		TFloat * const get_upbs(){
			return this->elkan_upper_base.get();
		}
		
		TFloat * const get_delta_C(){
			return this->elkan_delta_C.get();
		}
		
		std::function<void(TInt)> update_3v0_L_lowers_upper_where_changes_ati(){
			//TODO : neaten up and move out
			
			return [this](TInt ti){
				TInt data0 = (ti*this->ndata)/this->nthreads;
				TInt data1 = ((ti + 1)*this->ndata)/this->nthreads;

				TInt ndcalcs_local = 0;
				kmeans::sparse_update_L_lowers_upper_where_changes_3v0<TInt, TFloat>(this->ncentroids, this->dimension, data0, data1, *this->ptrdata, this->get_C(), this->get_data_l22s() + data0, this->get_C_l22s(), this->get_delta_C(), this->where_label_changes[ti], ndcalcs_local, this->get_L() + data0, this->get_lowers()  + data0*this->ncentroids, this->get_upbs() + data0);
				this->ndcalcs_X += ndcalcs_local;
			};
		}
		
		std::function<void(TInt)> update_S_H_from_where_changes_ati(){
			return [this](TInt ti){
					if (ti == 0){
						sparse::update_S_H_from_label_changes(*this->ptrdata, this->where_label_changes, this->get_sums(), this->get_counts(), this->nchanges);
					}
				};
		}
			
		
		
	public:
		typedef kmeans::BaseSparseElkan<TInt, TFloat> EB;
		template<typename... Args>
		SP3V0(Args&&... args): EB(std::forward<Args>(args)...)


		{
			this->setalgname("SP3V0");
			this->elkan_delta_C.reset(new TFloat [this->getncentroids()]);
		}
		virtual ~SP3V0(){}

		virtual TInt get_approximate_memory_requirement(){
			return EB::get_approximate_memory_requirement() + 
			sizeof(TFloat)*this->getncentroids(); // delta_C  
		}

		virtual void verbose_write_additional(){
			this->EB_verbose_write_additional();
			/* anything else to add ? */
		}

		virtual void set_initialisation_tasks(){
			/* all Elkan variants have same initialisation tasks */
			this->ElkBase_set_initialisation_tasks();
		}
	
		virtual void set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX)
			};
		}
		
		virtual void set_X_tasks(){
			this->X_tasks = {
				this->update_3v0_L_lowers_upper_where_changes_ati(),
				this->update_S_H_from_where_changes_ati()
			};
		}
};

}

#endif
