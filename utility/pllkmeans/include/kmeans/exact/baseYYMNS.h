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

#ifndef PLL_PLLYINYANGMNSBASEKMEANS_H
#define PLL_PLLYINYANGMNSBASEKMEANS_H

#include "baseYY.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{


template <typename TInt, typename TTau, typename TFloat>
/* " max  norm sum  " for lower bounds versions inherit from here */ 
class BaseYYMNS : public kmeans::BaseYY<TInt, TFloat>{
	
	private:
		std::unique_ptr<TTau []> tau;	

		std::unique_ptr<TTau []> tau_globallowers;		
		std::unique_ptr<TFloat []> globallowers_at_last; //ndata x 1
		
		TInt cutperiod;
		std::unique_ptr<TFloat []>  C_hist;
		std::unique_ptr<TFloat []>  C_l22s_hist; 
		std::unique_ptr<TFloat []>  u_delta_C;
		std::unique_ptr<TFloat []>  u_delta_G;
		//max change over all centroids:
		std::unique_ptr<TFloat []>  u_delta_global;

	protected:	
	
		TFloat * const get_globallowers_at_last(){
			return globallowers_at_last.get();
		}
		
		TFloat * const get_glowers_at_last(){
			return this->get_glowers_base();
		}
		
		TFloat * const get_upb_at_last(){
			return this->get_upb_base();
		}
		
		TTau * const get_tau(){
			return tau.get();
		}		
		
		TTau * const get_tau_globallowers(){
			return tau_globallowers.get();
		}		
		
		
		TFloat * get_C_hist(){
			return C_hist.get();
		}
		
		TFloat * get_C_l22s_hist(){
			return C_l22s_hist.get();
		}
		
		TFloat * const get_u_delta_C(){
			return u_delta_C.get();
		}
		
		TFloat * const get_u_delta_G(){
			return u_delta_G.get();
		}
		
		TFloat * const get_u_delta_global(){
			return u_delta_global.get();
		}
		

		
	/* with inspiration from pllelkankmeans4v2 */	
	public:
		typedef kmeans::BaseYY<TInt, TFloat> YYB;
		template<typename... Args>
		BaseYYMNS(Args&&... args): YYB(std::forward<Args>(args)...), 
		
		tau { new TTau [this->getndata()*this->get_ngroups()] }, 
		tau_globallowers { new TTau [this->getndata()*1] },
		
		globallowers_at_last { new TFloat [this->getndata()*1] },
				
		/* using memory derived cutperiod : P ncentroids dimension < max( ndata ngroups, ndata dimension ) --> 
		 * P < ndata/ncentroids max (ngroups / dimension, 1)
		 * (with an additional bounds 1 <= P <= 50 thrown on for good mesure)
		 * */
		cutperiod { std::min(static_cast<TInt>(50),
		std::max( static_cast<TInt>(1),
		static_cast<TInt>((static_cast<TFloat>(this->getndata())/static_cast<TFloat>(this->getncentroids())) * 
		std::max(static_cast<TFloat>(this->get_ngroups())/static_cast<TFloat>(this->getdimension()), static_cast<TFloat>(1))
		))) },
		C_hist { new TFloat [this->getncentroids()*this->getdimension()*(cutperiod)] },
		C_l22s_hist { new TFloat [this->getncentroids()*(cutperiod)] },
		u_delta_C { new TFloat [this->getncentroids()*(cutperiod + 1)] },
		u_delta_G { new TFloat [this->get_ngroups()*(cutperiod + 1)] },
		u_delta_global { new TFloat [1*(cutperiod + 1)] }
			
		{
			this->setalgname("BaseYYMNS");
			std::fill_n(tau.get(), this->getndata()*this->get_ngroups(), 0);
			std::fill_n(tau_globallowers.get(), this->getndata()*1, 0);
			std::fill_n(u_delta_C.get(), this->getncentroids()*(cutperiod + 1), 0);
			std::fill_n(u_delta_G.get(), this->get_ngroups()*(cutperiod + 1), 0);
			std::fill_n(u_delta_global.get(), 1*(cutperiod + 1), 0);
		}
		
		
		virtual TInt get_approximate_memory_requirement(){
			return YYB::get_approximate_memory_requirement() + 
			sizeof(TTau)*(
				this->getndata()*this->get_ngroups() + //tau
				this->getndata()) + //tau_globallowers
			sizeof(TFloat)*(
				this->getndata() + //globallowers_at_last
				this->getncentroids()*this->getdimension()*(cutperiod) + //C_hist
				this->getncentroids()*(cutperiod) + //C_l22s_hist
				this->getncentroids()*(cutperiod + 1) + //u_delta_C  
				this->get_ngroups()*(cutperiod + 1) + //u_delta_G
				(cutperiod + 1)); //u_delta_global
		}
				
			
		
		
		virtual ~BaseYYMNS(){}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to BaseYYMNS..\n\n";
		}

		
		TInt get_cutperiod(){
			return cutperiod;
		}
		
		void yinyang_MNS_initialisation_tasks(){
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_groups_glowers_upb_base_S_H_mati();
			/* add the task of setting globallowers_at_last. This c(sh)ould be merged the task where glowers is set, for 1 fewer barrier */
			this->initialisation_tasks.emplace_back(
				
				[this](TInt ti){
					TInt x0 = this->getndata()*ti/this->getnthreads();
					TInt x1 = this->getndata()*(ti + 1)/this->getnthreads();
					for (TInt i = x0; i < x1; ++i){
						globallowers_at_last[i] = *std::min_element(this->get_glowers_at_last() + i*this->get_ngroups(),  this->get_glowers_at_last() + (i+1)*(this->get_ngroups()));
					}
				}
				
			);		
		}		
		
		virtual void set_initialisation_tasks() = 0;
		
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
};

}

#endif

