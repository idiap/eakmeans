/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_ELKANKMEANS_MB3V0_H
#define PLL_ELKANKMEANS_MB3V0_H

#include "baseelkanminibatch.h"
#include "alg_X_selkSN.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class MB3V0 : public kmeans::BaseElkanMiniBatch<TInt, TFloat>{
	
	private: 
	
		//update L, lowers, upbs, S, H using data in range [x0, x1], perform on thread ti.
		void update_L_lowers_upper_S_H_3v0(TInt x0, TInt x1, TInt ti){
				
			if (this->round < this->mba.nsubrounds){
				this->set_upper_lowers_L(x0, x1);
				this->set_S_H(x0, x1, ti);
				this->mba.nchanges_on_batch[(this->mba.subround + 1)%this->mba.nsubrounds] = (x1 - x0);
			}
			
			else{
				this->mb_pll_principal_X(
				kmeans::update_L_lowers_upper_S_H_3v0<TInt, TFloat>, 
				ti,
				x1 - x0,
				this->data + x0*this->dimension,
				this->get_C(),
				this->get_data_l22s() + x0, 
				this->get_C_l22s(), 
				this->get_delta_C_fullround(), 
				this->get_L() + x0, 
				this->get_lowers() + x0*this->ncentroids, 
				this->get_upbs() + x0, 
				this->round);
			}
		}
	
		virtual std::function<void(TInt)> mb_update_L_lowers_upper_S_H_3v0_ati(){
			//TODO: this pattern should appear in baseminibatch.
			return [this](TInt ti){
				
				//the batch to use this round (same for all threads)
				TInt ind0 = this->mba.batchsize*((this->mba.subround + 1)%this->mba.nsubrounds); //as per simple;
				TInt ind1 = std::min(ind0 + this->mba.batchsize, this->ndata);
								
				//absolute indices of data to process on this threads
				TInt thisbatchsize = ind1 - ind0;
				TInt x0 = ind0 + (ti*thisbatchsize)/this->nthreads;
				TInt x1 = ind0 + ((ti + 1)*thisbatchsize)/this->nthreads;
			
				//TODO : implement the following
				//within: requires meta->delta_C of sorts. ndcalcs_X updated within
				//requires different behaviour on first view of minibatch.
				this->update_L_lowers_upper_S_H_3v0(x0, x1, ti); 	

			};
		}

			
	protected:
	
		std::unique_ptr<TFloat []> delta_C_fullround;
		std::unique_ptr<TFloat []> delta_C_history;
		
		TFloat * const get_lowers(){
			return this->elkan_lowers_base.get();
		}
		
		TFloat * const get_upbs(){
			return this->elkan_upper_base.get();
		}
		
		TFloat * const get_delta_C(){
			return this->elkan_delta_C.get();
		}
		
		TFloat * const get_delta_C_fullround(){
			return this->delta_C_fullround.get();
		}
		
		TFloat * const get_delta_C_history(TInt subround = 0){
			return this->delta_C_history.get() + subround*this->ncentroids;
		}
		
		virtual void verbose_write_additional(){
			
		}

		virtual void set_initialisation_tasks(){
			this->MBElkBase_set_initialisation_tasks(this->mba);
		}
	
		virtual void set_C_tasks(){

			this->C_tasks = {
				
				
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX),
				
				
				[this](TInt ti){
					if (ti == 0){

						arrutilv2::subtractfrom(this->ncentroids, this->get_delta_C_history(this->mba.subround), this->get_delta_C_fullround());
						arrutilv2::addto(this->ncentroids, this->get_delta_C(), this->get_delta_C_fullround());

						//update delta_C_subrounds (move delta_C into delta_C_subrounds)
						std::memcpy(this->get_delta_C_history(this->mba.subround), this->get_delta_C(), sizeof(TFloat)*this->ncentroids);
					}
				}
			};			
		}
		
		virtual void set_X_tasks(){
			this->X_tasks = {
				//[this](TInt ti){
					//std::cout << "\nnchanges_on_batch" << std::endl;
					//for (TInt a = 0; a < this->mba.nsubrounds; ++a){
						//std::cout << this->mba.nchanges_on_batch[a] << " ";
					//}
					//std::cout << std::endl;
				//},
				
				this->mb_update_L_lowers_upper_S_H_3v0_ati(),
				this->minibatch_subround_update(this->mba)
			};
		}
		
		
	public:
		typedef kmeans::BaseElkanMiniBatch<TInt, TFloat> EB;
		template<typename... Args>
		MB3V0(Args&&... args): EB(std::forward<Args>(args)...)


		{
			this->setalgname("MB3V0");
			this->elkan_delta_C.reset(new TFloat [this->getncentroids()]);
			this->delta_C_fullround.reset(new TFloat [this->ncentroids]);
			this->delta_C_history.reset(new TFloat [this->ncentroids * this->mba.nsubrounds]);
			std::fill_n(this->get_delta_C_fullround(), this->ncentroids, 0);
			std::fill_n(this->get_delta_C_history(), this->mba.nsubrounds*this->ncentroids, 0);			
		}
		virtual ~MB3V0(){}
		
		virtual TInt get_approximate_memory_requirement(){
			return EB::get_approximate_memory_requirement() + 
			sizeof(TFloat)*this->getncentroids()*(2 + this->mba.nsubrounds); // delta_C etc   
		}
};

}

#endif


