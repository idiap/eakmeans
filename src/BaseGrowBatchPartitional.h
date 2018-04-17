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

#ifndef BASEGROWBATCHPARTITIONAL_H
#define BASEGROWBATCHPARTITIONAL_H

#include "BaseGrowBatch.h"

namespace kmeans{
template <typename TInt, typename TFloat>
/* A type of GrowBatch, so batch size doubles when determined to be appropriate. Specifically, doubling is performed when the l2 norm between center (vector R^{dk}) in partition A and B (random partitions of data) are larger than the l2 norm of current center to center in previous round  */
class BaseGrowBatchPartitional : public kmeans::BaseGrowBatch<TInt, TFloat>{
	
	private:		

								
		void BGBP_constructor_helper(){
			
			this->setalgname("Base Grow Batch Partitional");
					
			this->C_A.reset(new TFloat [this->ncentroids*this->dimension] );
			this->C_A_l22s.reset(new TFloat [this->ncentroids] );	
			this->sums_A.reset(new TFloat [this->ncentroids*this->dimension] );
			this->counts_A.reset(new TInt [this->ncentroids] );
			
			this->C_B.reset(new TFloat [this->ncentroids*this->dimension] );
			this->C_B_l22s.reset(new TFloat [this->ncentroids] );	
			this->sums_B.reset(new TFloat [this->ncentroids*this->dimension] );
			this->counts_B.reset(new TInt [this->ncentroids] );

			std::fill_n(this->sums_A.get(), this->ncentroids*this->dimension, 0);
			std::fill_n(this->counts_A.get(), this->ncentroids, 0);
			std::fill_n(this->sums_B.get(), this->ncentroids*this->dimension, 0);
			std::fill_n(this->counts_B.get(), this->ncentroids, 0);

			this->delta_C2_AB.reset(new TFloat [this->ncentroids]);
			
		}
	
	protected:
			
		virtual void set_initialisation_tasks() = 0 ;

		virtual bool should_double() override final{
			return (this->delta_C2_AB_scalar > this->gba.threshold*delta_C2_scalar);
		}


		//ever used? Looks buggy (this->ndata??)
		template <typename Function, typename... Args>
		void pll_principal_X(const Function & X_updater, TInt ti, const char & partition, Args&&... args){
			TInt data_start = this->get_partition_offset(partition);
			TInt data_end = data_start + this->gba.ndata_active/2;
			TInt x0 = data_start + (ti*this->ndata)/this->nthreads;
			TInt x1 = data_start + ((ti+1)*this->ndata)/this->nthreads;		

			arrutilv2::pll_update_L_etc(X_updater, 
			this->ncentroids, this->dimension, this->get_partition_sums(partition), this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_partition_counts(partition), this->get_dcounts() + ti*this->ncentroids, this->nchanges, this->ndcalcs_X, this->work_mutex, x1-x0, this->data +x0*this->dimension, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), 
			std::forward<Args>(args)...);

		}


		std::function<void(TInt)> GBP_set_S_H_ati_on_partition(const char & partition){
			return [partition, this](TInt ti){
				TInt ndatatouse = this->gba.ndata_active/2;
				TInt x0 = get_partition_offset(partition) + (ti*ndatatouse)/this->nthreads;
				TInt x1 = get_partition_offset(partition) + ((ti+1)*ndatatouse)/this->nthreads;
	
				arrutilv2::set_S_H(x1-x0, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_partition_sums(partition), this->get_partition_counts(partition), this->work_mutex);
			};
		}
		
		std::vector<std::function<void(TInt)>> GBP_set_S_H_mati(){
			
			std::vector<std::function<void(TInt)>> tasks;
			//set partitional sums, counts
			tasks.push_back(GBP_set_S_H_ati_on_partition('A'));
			tasks.push_back(GBP_set_S_H_ati_on_partition('B'));
			
			
			//set global sums, counts
			tasks.push_back(arrutilv2::set_vector_sum_ati(this->nthreads, this->dimension*this->ncentroids, this->get_partition_sums('A'), this->get_partition_sums('B'), this->get_sums()));

			tasks.push_back(arrutilv2::set_vector_sum_ati(this->nthreads, this->ncentroids, this->get_partition_counts('A'), this->get_partition_counts('B'), this->get_counts()));
		
			return tasks;
		
		}
		
			
	
		
		/* centroid specifics for first (A) and second (B) halves of ndata_active */
		std::unique_ptr<TFloat []> C_A;
		std::unique_ptr<TFloat []> C_A_l22s;
		std::unique_ptr<TFloat []> sums_A;
		std::unique_ptr<TInt []> counts_A; 
	 
		std::unique_ptr<TFloat []> C_B;
		std::unique_ptr<TFloat []> C_B_l22s;
		std::unique_ptr<TFloat []> sums_B;
		std::unique_ptr<TInt []> counts_B; 
				
		/* used to determine if exapansion should take place (and maybe other things) */
		std::unique_ptr<TFloat []> delta_C2_AB;
		
		/*  \|C_{t} - C_{t-1}\|_2  */
		TFloat delta_C2_scalar; 
		
		/*  \|C_A - C_B\|_2  */
		TFloat delta_C2_AB_scalar; 


		// MUST set delta_C in this function because set_C_tasks does not do it itself.
		virtual std::vector<std::function<void(TInt)>> specific_C_tasks_mati() = 0;

		// This function will probably have an 'A' and 'B' X-update specific function.
		virtual std::vector<std::function<void(TInt)>> specific_X_tasks_mati() = 0;


		inline TInt get_partition_offset(const char & partition){
			if (partition == 'A'){
				return 0;
			}
			else if (partition == 'B'){
				return this->gba.ndata_active/2;
			}
			else {
				throw std::logic_error("partition in get_partition_offset should be 'A' or 'B', not" + std::to_string(partition));
			}
		}


		inline TFloat * const get_partition_sums(const char & partition){
			
			if (partition == 'A'){
				return this->sums_A.get();
			}
			
			else if (partition == 'B'){
				return this->sums_B.get();
			}
			
			else {
				throw std::logic_error("partition in get_partition_sums should be 'A' or 'B', not" + std::to_string(partition));
			}
			
		}

		inline TFloat * const get_partition_C(const char & partition){
			
			if (partition == 'A'){
				return this->C_A.get();
			}
			
			else if (partition == 'B'){
				return this->C_B.get();
			}
			
			else {
				throw std::logic_error("partition in get_partition_C should be 'A' or 'B', not" + std::to_string(partition));
			}
			
		}
		
		//TODO : do I use this everywhere I should?
		inline TFloat * const get_partition_C_l22s(const char & partition){
			
			if (partition == 'A'){
				return this->C_A_l22s.get();
			}
			
			else if (partition == 'B'){
				return this->C_B_l22s.get();
			}
			
			else {
				throw std::logic_error("partition in get_partition_C_l22s should be 'A' or 'B', not" + std::to_string(partition));
			}
			
		}
		
		
		
		inline TInt * const get_partition_counts(const char & partition){
			
			if (partition == 'A'){
				return this->counts_A.get();
			}
			
			else if (partition == 'B'){
				return this->counts_B.get();
			}
			
			else {
				throw std::logic_error("partition in get_partition_counts should be 'A' or 'B', not" + std::to_string(partition));
			}
		}
		
		
		
		inline TFloat * const get_delta_C2_AB(){
			return this->delta_C2_AB.get();
		}






		virtual void set_X_tasks() final{
			
			//If just had an increase in ndata_active, combine sums, counts A, B. TODO : move from GBSimple to here. 
			
			//Do the basic 'A' and 'B' updates, partition specific sums and counts are updated.
			this->X_tasks = specific_X_tasks_mati();
			
			//Update global sums and counts
			this->X_tasks.push_back(arrutilv2::set_vector_sum_ati(this->nthreads, this->dimension*this->ncentroids, this->get_partition_sums('A'), this->get_partition_sums('B'), this->get_sums()));

			this->X_tasks.push_back(arrutilv2::set_vector_sum_ati(this->nthreads, this->ncentroids, this->get_partition_counts('A'), this->get_partition_counts('B'), this->get_counts()));
			
			//this->X_tasks.push_back([this](TInt ti){if (ti == 0) { std::cout << "`end of X  tasks " << std::endl; } } );


		}
		

		virtual void set_C_tasks() final{
			
			//class specific tasks (should update C, C_l22s, delta_C, etc.)
			this->C_tasks = specific_C_tasks_mati();
			
			//the usual GB suspects,
			this->C_tasks.push_back(

				//set C_A C_A_l22s
				arrutilv2::update_C_C_l22s_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_partition_sums('A'), this->get_partition_counts('A'), this->get_partition_C('A'), this->get_partition_C_l22s('A')));
				
			
			this->C_tasks.push_back(
			
				//set C_B C_B_l22s
				arrutilv2::update_C_C_l22s_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_partition_sums('B'), this->get_partition_counts('B'), this->get_partition_C('B'), this->get_partition_C_l22s('B')));

			
			this->C_tasks.push_back(
				//set delta_C2_AB (vector)
				arrutilv2::set_rl22s_ati(this->nthreads, this->ncentroids, this->dimension,  
				this->get_partition_C('A'), 
				this->get_partition_C('B'), 
				this->get_delta_C2_AB(), 
				this->ndcalcs_X));
			
			this->C_tasks.push_back(
				//save d_C__over__d_AB (hacky for output) and 
				//set delta_C2_scalar and delta_C2_AB
				[this](TInt ti){
					if (ti == 0){
						//std::cout << ">>>>>>>" << this->delta_C2_scalar/std::max<TFloat>(1e-7, this->delta_C2_AB_scalar) << "\t";
						delta_C2_scalar = 0;
						delta_C2_AB_scalar = 0;
						
						for (TInt ci = 0; ci < this->ncentroids; ++ci){
							delta_C2_scalar += this->gba.delta_C[ci]*this->gba.delta_C[ci];
							delta_C2_AB_scalar += this->delta_C2_AB[ci];
						}
						
						//I don't think this is used, see basecluster.h						
						delta_C2_scalar = std::sqrt(std::max(static_cast<TFloat>(0.), delta_C2_scalar));
						delta_C2_AB_scalar = std::sqrt(std::max(static_cast<TFloat>(0.), delta_C2_AB_scalar));
						this->gba.d_C__over__d_AB = this->delta_C2_scalar/std::max<TFloat>(1e-7, this->delta_C2_AB_scalar);

					}
				}
			);
				

			this->C_tasks.push_back(
				//store previous size. test, check if size should increase.				
				this->update_ndata_active_ati(this->gba)
			);
			
				//TODO: what to do with sums and counts in X-step if increase? I suppose partition A takes all, B from scratch (Unless less than double)  
			
	 }

	public:

		template<typename... Args>
		BaseGrowBatchPartitional(TInt batchsize0, Args&&... args): kmeans::BaseGrowBatch<TInt, TFloat> (batchsize0, std::forward<Args>(args)...) //1, 21		
		{
			this->BGBP_constructor_helper();
			throw std::runtime_error("It is possible that the change I made in BGB_constructor_helper is going to break this algorithm : changed gbapp.ndata_active = 2*batchsize0; to gbapp.ndata_active = batchsize0;.");
		}
			
		virtual ~BaseGrowBatchPartitional(){};
	

};
}



#endif
