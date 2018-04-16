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

#ifndef BASEGROWBATCH_H
#define BASEGROWBATCH_H

#include "basekmeans.h"

#include "growbatchapp.h"


namespace kmeans{
template <typename TInt, typename TFloat>
class BaseGrowBatch : public kmeans::BaseKmeans<TInt, TFloat>{

	private:
		virtual void endroundupdate() override final{
			this->iscomplete = (this->nchanges == 0 && (this->gba.ndata_active == this->ndata)) || (this->duration > this->maxtime) || (this->round >= this->maxrounds);
			this->nchanges = 0;
			++this->round;
		}

	protected:	
	
		growbatchapp::GBApp<TInt, TFloat> gba;
		//For the data used in an X update, some algorithms 
		//do full on data x centroid multiplications, how much 
		//data can be used per thread in a full data x centroid
		// product?
		TInt maxpermultiplyblock;	

		virtual void set_X_tasks() = 0;
		virtual void set_C_tasks() = 0;
		
		inline TFloat * const get_delta_C(){
			return this->gba.delta_C.get();
		}
		
		virtual void set_initialisation_tasks() = 0;		

		void BGB_constructor_helper_densebits(){

			//stuff specific to dense goes here.
			
			this->setalgname("Dense Base Grow Batch");
			
			this->maxpermultiplyblock =
			std::max(static_cast<TInt> (1),
			static_cast<TInt> ((this->getndata() * this->getdimension())/(2 * this->getncentroids() * this->nthreads)));			
		}
		
		
		//A c&p from minibatch base. Not as code reducing as the baseexact version, but easier to understand
		template <typename Function, typename... Args>
		void gb_pll_principal_X(const Function & X_updater, TInt ti, Args&&... args){
	
			arrutilv2::pll_update_L_etc(
			//The compulsory parameters to pll_update_L_etc,
			X_updater, 
			this->ncentroids, this->dimension, this->get_sums(), this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_counts(), this->get_dcounts() + ti*this->ncentroids, this->nchanges, this->ndcalcs_X, this->work_mutex
			//The additional parameters to pll_update_L_etc with correct offset
			, std::forward<Args>(args)...);
		}
		
		
	private:
		virtual void set_summaries() override final{
			this->set_summaries_growbatch(this->gba);
		}


		//Note that BaseGrowBatchMse overrides this.
		virtual void set_mse() override {
			//if not all data is active, refuse to compute the mse
			if (this->gba.ndata_active != this->ndata){
				this->mse = -1;
			}
			
			else{
				this->mse = arrutilv2::getmeanl22at(this->ncentroids, this->dimension, this->get_C(), this->ndata, this->data, this->get_L(), this->get_C_l22s(), this->get_data_l22s());
			}
		}
		
	
	
		
		
		

	public:
	
		template<typename... Args>
		BaseGrowBatch(TInt batchsize0, Args&&... args): kmeans::BaseKmeans<TInt, TFloat> (std::forward<Args>(args)...)		
		{
			this->BGB_constructor_helper(batchsize0, this->gba);
			this->BGB_constructor_helper_densebits();
		}
		virtual ~BaseGrowBatch(){};
};
}

#endif
