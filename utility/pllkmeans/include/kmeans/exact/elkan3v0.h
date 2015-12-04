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

#ifndef PLL_ELKANKMEANS_3V0_H
#define PLL_ELKANKMEANS_3V0_H

#include "baseelkan.h"
#include "alg_X_selkSN.h"

namespace kmeans{

/* discrepency in ndcalcs as compared to a3v0 due to not computing CC initially (I propose) */

template <typename TInt, typename TFloat>
class P3V0 : public kmeans::BaseElkan<TInt, TFloat>{
	
	private:
		std::unique_ptr<TFloat []> delta_C;
		
	protected:
		TFloat * const get_lowers(){
			return this->get_lowers_base();
		}
		
		TFloat * const get_upbs(){
			return this->get_upper_base();
		}
		
		TFloat * const get_delta_C(){
			return delta_C.get();
		}
		
		std::function<void(TInt)> update_L_lowers_upbs_S_H_3v0_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				this->pll_principal_X(update_L_lowers_upbs_S_H_3v0<TInt, TFloat>, ti, this->get_delta_C(), this->get_L() + x0,  this->get_lowers() + x0*this->getncentroids(), this->get_upbs() + x0, this->round);
			};
		}
		
		
	public:
		typedef kmeans::BaseElkan<TInt, TFloat> EB;
		template<typename... Args>
		P3V0(Args&&... args): EB(std::forward<Args>(args)...), 

		delta_C{ new TFloat [this->getncentroids()] }
		{
			this->setalgname("p3v0");
			//std::cout << "memory request " << this->get_approximate_memory_requirement() / (1024.*1024. *1024.) << " GB " << std::endl;
		}
		virtual ~P3V0(){}

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
			this->EB_set_initialisation_tasks();
		}
	
		virtual void set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX)
			};
		}
		
		virtual void set_X_tasks(){		
			this->X_tasks = {
				this -> update_L_lowers_upbs_S_H_3v0_ati()
			};
		}
};

}

#endif
