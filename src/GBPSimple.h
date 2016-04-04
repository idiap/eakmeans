/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef GBPSimple_H
#define GBPSimple_H

#include "BaseGrowBatchPartitional.h"

namespace kmeans{
template <typename TInt, typename TFloat>
class GBPSimple : public kmeans::BaseGrowBatchPartitional<TInt, TFloat>{
	
	private: 
	
		/* I don't think this needs to be split by partition */
		std::function<void(TInt)> set_L_ati(const char & partition){
		
			return [this, partition](TInt ti){
				TInt local_ndcalcs = 0;
				TInt x0 = this->get_partition_offset(partition) + (ti*this->gba.ndata_active/2)/this->nthreads;
				TInt x1 = this->get_partition_offset(partition) + ((ti+1)*this->gba.ndata_active/2)/this->nthreads;				
				arrutilv2::set_rargmins(x1 - x0, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, local_ndcalcs);				
				this->ndcalcs_X += local_ndcalcs;
				this->nchanges += x1 - x0;
			};
		}
			
	protected:


		
		inline void update_L_S_H_batch(const char & partition, TInt ti, TInt x0, TInt x1){
			
			arrutilv2::update_L_S_H_batch(x1-x0, this->maxpermultiplyblock, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_partition_sums(partition), this->get_partition_counts(partition), this->nchanges, this->work_mutex);
		}
		
		
		inline void set_L_S_H_batch(const char & partition, TInt ti, TInt x0, TInt x1){
			arrutilv2::update_L_S_H_batch_increment_only(x1-x0, this->maxpermultiplyblock, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_partition_sums(partition), this->get_partition_counts(partition), this->nchanges, this->work_mutex);
		}

		
		
		std::function<void(TInt)> prepare_for_partition_update_L_S_H_batch_ati(){
			
			return [this](TInt ti){
				if (ti == 0){
	
					//size did not increase, nothing to do
					if (this->gba.ndata_active == this->gba.ndata_active_previous){
							
					}
					
					else{
						//size exactly doubled. 
						if (this->gba.ndata_active == 2*this->gba.ndata_active_previous){
							//'A takes A + B'
							std::memcpy(this->get_partition_sums('A'), this->get_sums(), sizeof(TFloat)*this->dimension*this->ncentroids);
							std::memcpy(this->get_partition_counts('A'), this->get_counts(), sizeof(TInt)*this->ncentroids);
							
							//B gets set to zero
							std::fill_n(this->get_partition_sums('B'), this->dimension*this->ncentroids, 0);
							std::fill_n(this->get_partition_counts('B'), this->ncentroids, 0);
								
						}
						
						//size less than doubled, so at capacity now. 
						else{
							//A and B get set to zero
							std::fill_n(this->get_partition_sums('A'), this->dimension*this->ncentroids, 0);
							std::fill_n(this->get_partition_counts('A'), this->ncentroids, 0);
							std::fill_n(this->get_partition_sums('B'), this->dimension*this->ncentroids, 0);
							std::fill_n(this->get_partition_counts('B'), this->ncentroids, 0);						
						}
					}	
				}			
			};	
		}
		
		
				





		
		inline std::function<void(TInt)> partition_update_L_S_H_batch_ati(const char & partition){
			
			return [this, partition](TInt ti){
				TInt x0 = this->get_partition_offset(partition) + (ti*this->gba.ndata_active/2)/this->nthreads;
				TInt x1 = this->get_partition_offset(partition) + ((ti+1)*this->gba.ndata_active/2)/this->nthreads;	

				//size did not increase, standard update ensues
				if (this->gba.ndata_active == this->gba.ndata_active_previous){
					update_L_S_H_batch(partition, ti, x0, x1);					
				}
				
				
				else{
					//size exactly doubled. 
					if (this->gba.ndata_active == 2*this->gba.ndata_active_previous){
						if (partition == 'A'){
							update_L_S_H_batch(partition, ti, x0, x1);
						}
						else{
							set_L_S_H_batch(partition, ti, x0, x1);
						}
						
					}
					
					//size less than doubled, so at capacity now. 
					else{
						set_L_S_H_batch(partition, ti, x0, x1);
					}
				}
				this->ndcalcs_X += (x1 - x0)*this->ncentroids;

				
			};
			
		}
		
		
		
		
		virtual std::vector<std::function<void(TInt)>> specific_X_tasks_mati() override{
	
			return {

				this->prepare_for_partition_update_L_S_H_batch_ati(),
				this->partition_update_L_S_H_batch_ati('A'),
				this->partition_update_L_S_H_batch_ati('B')
			};					
			
		}
		
		
		virtual std::vector<std::function<void(TInt)>> specific_C_tasks_mati(){
			std::vector<std::function<void(TInt)>> specific_C_tasks;

			specific_C_tasks.push_back(			
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->gba.delta_C.get(), this->ndcalcs_notX)
			);
			return specific_C_tasks;
		}
				
		//some initialisation, using the first batch if necessary
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_inds0_mati(){
			return this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), this->gba.ndata_active);
		}
		
		//use makeset_C_C_l22s_inds0_mati and set L of first batch (all others to 0)
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_mati(){
			
			std::vector<std::function<void(TInt)> > tasks;
			tasks = this->makeset_C_C_l22s_inds0_mati();			
			
			//set all L to "ncentroids" so that it is reflected that all change in zeroth round
			tasks.push_back(
				[this](TInt ti){
					if (ti == 0){
						std::fill_n(this->get_L(), this->ndata, this->ncentroids);
					}
				}
			);
			
			//set labels on active data
			tasks.push_back(
				this->set_L_ati('A')
			);
			
			tasks.push_back(
				this->set_L_ati('B')
			);
				
			return tasks;
		}
		
		virtual void set_initialisation_tasks(){
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_mati();
			auto init_tasks_B = this->GBP_set_S_H_mati();	
			this->initialisation_tasks.insert(this->initialisation_tasks.end(), init_tasks_A.begin(), init_tasks_A.end());
			this->initialisation_tasks.insert(this->initialisation_tasks.end(), init_tasks_B.begin(), init_tasks_B.end());

		}

	

	public:
		template<typename... Args>
		GBPSimple(Args&&... args): kmeans::BaseGrowBatchPartitional<TInt, TFloat> (std::forward<Args>(args)...)	
		{

		}
		
		
		virtual ~GBPSimple(){};
};
}

//extern template class kmeans::GBPSimple<size_t, double>;
//extern template class kmeans::GBPSimple<size_t, float>;



#endif

