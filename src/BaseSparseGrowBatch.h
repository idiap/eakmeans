/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

//Based blindly on BaseSparseGrowBatch.h

#ifndef BASESPARSEGROWBATCH_H
#define BASESPARSEGROWBATCH_H

#include "basesparsekmeans.h"
#include "growbatchapp.h"


namespace kmeans{
template <typename TInt, typename TFloat>
class BaseSparseGrowBatch : public kmeans::BaseSparseKmeans<TInt, TFloat>{

	private: 	
	
		virtual void endroundupdate() override final{
			this->iscomplete = (this->nchanges == 0 && (this->gba.ndata_active == this->ndata)) || (this->duration > this->maxtime) || (this->round >= this->maxrounds);
			this->nchanges = 0;
			++this->round;
		}
	
		virtual void sgb_update_L_etc(TInt x0, TInt x1, TInt ti) = 0;
		virtual void sgb_set_L_etc(TInt x0, TInt x1, TInt ti) = 0;

	protected:	
	
		growbatchapp::GBApp<TInt, TFloat> gba;

		virtual void set_C_tasks() = 0;
		
		inline TFloat * const get_delta_C(){
			return this->gba.delta_C.get();
		}
		
		virtual void set_initialisation_tasks() = 0;		

		void BGB_constructor_helper_sparsebits(){
			this->setalgname("Sparse Base Grow Batch");
		}
		
		
		//This may be a bit premature : hoping that my assumption of one task per update round is ~accurate.
		virtual void set_X_tasks() override final{ 
			this->X_tasks = this->sgb_update_L_dn_etc_S_H_mati(); 		
		}

		
		std::vector<std::function<void(TInt)> > sgb_update_L_dn_etc_S_H_mati(){
		std::vector<std::function<void(TInt)> > tasks = {};
			tasks.emplace_back (
				//update L and dn of data used in previous round			
				[this](TInt ti){
					TInt x0 = (ti*this->gba.ndata_active_previous)/this->nthreads;
					TInt x1 = ((ti+1)*this->gba.ndata_active_previous)/this->nthreads;
					
					this->sgb_update_L_etc(x0, x1, ti);
					
				}
			);
			
			tasks.emplace_back(
				[this](TInt ti){
					if (ti == 0){
						sparse::update_S_H_from_label_changes(*this->ptrdata, this->where_label_changes, this->get_sums(), this->get_counts(), this->nchanges);	 
					}
				}
			);
				
			tasks.emplace_back(
				[this](TInt ti){
				//set L and dn of unused 
					if (this->gba.ndata_active != this->gba.ndata_active_previous){
						TInt ndata_tail = this->gba.ndata_active - this->gba.ndata_active_previous;
						TInt x0 =  this->gba.ndata_active_previous + (ti*ndata_tail)/this->nthreads;
						TInt x1 =  this->gba.ndata_active_previous + ((ti + 1)*ndata_tail)/this->nthreads;
						this->sgb_set_L_etc(x0, x1, ti);
						this->nchanges += x1 - x0;
					} 
				}
			);
			
			tasks.emplace_back(
				[this](TInt ti){
					if (ti == 0){
						if (this->gba.ndata_active != this->gba.ndata_active_previous){
							sparse::increment_S_H(this->gba.ndata_active_previous, 
							this->gba.ndata_active, 
							*this->ptrdata, 
							this->get_L(), 
							this->get_sums(), 
							this->get_counts());
						}
					}
				}
			);
			return tasks;
		}
		
		
	private:
		virtual void set_summaries() override final{
			this->set_summaries_growbatch(this->gba);
		}

		virtual void set_mse() override {
			if (this->gba.ndata_active != this->ndata){
				this->mse = -1;
			}
			
			else{
				this->mse = this->getmeanl22at();
			}
		}
		
	
	
		
		
		

	public:
	
		template<typename... Args>
		BaseSparseGrowBatch(TInt batchsize0, Args&&... args): kmeans::BaseSparseKmeans<TInt, TFloat> (std::forward<Args>(args)...)
		{
			this->BGB_constructor_helper(batchsize0, this->gba);
			this->BGB_constructor_helper_sparsebits();

		}
		virtual ~BaseSparseGrowBatch(){};
};
}

#endif
