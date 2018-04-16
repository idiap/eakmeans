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

#ifndef PLL_BASEDENSECENTROIDKMEANSSS_H
#define PLL_BASEDENSECENTROIDKMEANSSS_H

#include "basecluster.h"
#include <cstring>
#include "arrutilv2l1.h"
#include "arrutilv2copy.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseDenseCentroidKmeans : public cluster::BaseCluster<TInt, TFloat> { 

	
	private:
		virtual void set_mse() = 0;
		virtual void set_summaries() = 0;




	protected:

		
		/* nothing to do with some classes, here cos it is used by more than 1 inherited class. This is a vile hack.*/
		//base elkan
		std::unique_ptr<TFloat []> elkan_lowers_base;
		std::unique_ptr<TFloat []> elkan_upper_base;	
		
		//elkan 3v0
		std::unique_ptr<TFloat []> elkan_delta_C;

		/*  *   *  * * ** *  *   *    */


		const TFloat * C_init;
		std::unique_ptr<TFloat []> data_l22s;
		std::unique_ptr<TFloat []> C;
		std::unique_ptr<TFloat []> C_l22s;
		std::unique_ptr<TFloat []> sums;
		std::unique_ptr<TInt []> counts; 
		
		std::unique_ptr<TFloat []> valdata_l22s;
	
		virtual std::unique_ptr<TFloat []> getnew_data_l22s() = 0;
		virtual std::unique_ptr<TFloat []> getnew_valdata_l22s() = 0;
	
		virtual void set_L(TInt x0, TInt x1) = 0;
		virtual void set_L_dn(TInt x0, TInt x1){
			throw std::logic_error("Call to set_L_dn (in basedensecentroidkmeans) is unexpected");
		}
		
		//used by Elkan kmeans
		virtual void set_upper_lowers_L(TInt x0, TInt x1){
		 throw std::logic_error("only overriden virtual versions of set_upper_lowers_L (in basedensecentroidkmeans.h) expected to be called, bailing. ");	
		};
		
		virtual void set_L_lowers_dn(TInt x0, TInt x1){
		 throw std::logic_error("only overriden virtual versions of set_L_uppers_dn (in basedensecentroidkmeans.h) expected to be called, bailing. ");			
		}

		std::function<void(TInt)> set_L_ati(TInt data0, TInt data1){
			return [this, data0, data1](TInt ti){

				TInt local_ndata = data1 - data0;
				TInt x0 = (data0 + ti*local_ndata)/this->nthreads;
				TInt x1 = (data0 + (ti+1)*local_ndata)/this->nthreads;
				
				this->set_L(x0, x1);
			};
		}
		
		std::function<void(TInt)> set_L_dn_ati(TInt data0, TInt data1){
			return [this, data0, data1](TInt ti){

				TInt local_ndata = data1 - data0;
				TInt x0 = (data0 + ti*local_ndata)/this->nthreads;
				TInt x1 = (data0 + (ti+1)*local_ndata)/this->nthreads;
				
				this->set_L_dn(x0, x1);
			};
		}

		
		std::function<void(TInt)> set_L_lowers_dn_ati(TInt data0, TInt data1){
			return [this, data0, data1](TInt ti){

				TInt local_ndata = data1 - data0;
				TInt x0 = (data0 + ti*local_ndata)/this->nthreads;
				TInt x1 = (data0 + (ti+1)*local_ndata)/this->nthreads;
				
				this->set_L_lowers_dn(x0, x1);
			};
		}

		virtual std::function<void(TInt)> base_set_S_H_ati(TInt data0, TInt data1) = 0;
		
		
		void assignmemory_elkan_upper_lowers(){
			this->elkan_lowers_base.reset(new TFloat [this->getndata()*this->getncentroids()]);
			this->elkan_upper_base.reset(new TFloat [this->getndata()]);
		}
		
		//useful for grow batch algorithms where uppers is taken care of as already have exact distances. 
		void assignmemory_elkan_lowers(){
			this->elkan_lowers_base.reset(new TFloat [this->getndata()*this->getncentroids()]);
		}
		
		
		
		TInt get_elkan_base_memory(){
			return
			sizeof(TFloat)*this->getndata()*this->getncentroids() + //lower base
			sizeof(TFloat)*this->getndata(); //upper base
		}
		
		

		//to be used by certain exact kmeans
		std::vector<std::function<void(TInt)> > exact_makeset_C_C_l22s_inds0_mati(){
			return this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), this->ndata);
		}


		//to be used by minibatchelkans
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati(TInt data0, TInt data1){
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_lowers_upper_mati(data0, data1);
			auto init_task_B = this->base_set_S_H_ati(data0, data1);
			auto initialisation_tasks = std::move(init_tasks_A);
			initialisation_tasks.push_back(std::move(init_task_B));	
			return initialisation_tasks;			
		}


		//to be used by exact elkans 
		std::vector<std::function<void(TInt)> > exact_makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati(){
			auto init_tasks_A = this->exact_makeset_C_C_l22s_L_inds0_lowers_upper_mati();
			auto init_task_B = this->base_set_S_H_ati(static_cast<TInt>(0), this->ndata);
			
			
			auto initialisation_tasks = std::move(init_tasks_A);
			initialisation_tasks.push_back(std::move(init_task_B));	
			return initialisation_tasks;
		}
		
		
		void ElkBase_set_initialisation_tasks(){
			this->initialisation_tasks = this->exact_makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati();
		}
		
		
		
		void MBElkBase_set_initialisation_tasks(minibatchapp::MiniBatchApp<TInt> & mba_in){
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati(static_cast<TInt>(0), mba_in.initialising_batch_size);
			//this->nchanges += mba_in.initialising_batch_size*this->ncentroids;
			
			
			//see minibatch_makeset_C_C_l22s_L_inds0_mati
			this->initialisation_tasks.push_back(
			[this, &mba_in](TInt ti){
				if (mba_in.nsubrounds > 1){
					mba_in.nchanges_on_batch[0] += mba_in.initialising_batch_size;
				}
			});
		}


		
		//{
				//return [](TInt){};
			//}// = 0;
		
		//TODO: I'd prefer not to have any functions over complete range here.
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_mati(){
			//will not be used by all who inherit	 (only simplest and basesimple)
			std::vector<std::function<void(TInt)> > tasks;
			tasks = this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), this->ndata);
			tasks.push_back(
				set_L_ati(static_cast<TInt> (0), this->ndata)
			);
			return tasks;
		}


				
		virtual void set_all_initialisation_tasks() final{
			
			this->all_initialisation_tasks = {
				
				[this](TInt ti){
					if (ti == 0){
						//TODO this is not optimal, should simply move in, but as many functions are initialised frozen to address, I do not want to change address of data_l22s. Embarrasing code (Is preceding true?)
						
						auto new_data_l22s = this->getnew_data_l22s();
						std::memcpy(this->data_l22s.get(), new_data_l22s.get(), sizeof(TFloat)*this->ndata);
						
						if (this->nvaldata > 0){
							auto new_data_l22s = this->getnew_valdata_l22s();
							std::memcpy(this->valdata_l22s.get(), new_data_l22s.get(), sizeof(TFloat)*this->nvaldata); //memcpy is SOOO dangerous! 2 hrs wasted because I had ndata instead of nvaldata here :() :( ) :(  )
				
						}
						
					}
				}
			
			};
			
			this->set_initialisation_tasks();
			this->all_initialisation_tasks.insert(this->all_initialisation_tasks.end(), this->initialisation_tasks.begin(), this->initialisation_tasks.end());
			
			//this->all_initialisation_tasks.push_back([](TInt ti){std::cout << "all initialisation tasks complete" << std::endl;});
			
			
		}
		
		virtual void set_initialisation_tasks() = 0;

		virtual std::function<void(TInt)> set_C_Cl22s_from_inds0etc_ati() = 0;		
		virtual std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt [] > > get_C_inds0_uniform(TInt ind0, TInt ind1) = 0;		
		virtual void do_kmeanspp_initialisation(TInt ind0, TInt ind1) = 0;

	
	
		std::vector<std::function<void(TInt)> > base_makeset_C_C_l22s_inds0_mati(TInt ind0, TInt ind1){ 
			
			std::vector<std::function<void(TInt)> > tasks;

			tasks.push_back(
					[this](TInt ti){
						if (ti == 0){

						}
					}
				);
				
			
			if (this->C_init != nullptr && this->initialisation_method.compare("from_C") == 0){
				this->inds0.reset(nullptr);				
				std::memcpy(this->get_C(), this->C_init, this->ncentroids*this->dimension*sizeof(TFloat));
				arrutilv2::set_rl22s(this->ncentroids, this->dimension, this->get_C(), this->get_C_l22s());
			}
			
			else if (this->data_indices_init_from != nullptr && this->initialisation_method.compare("from_indices") == 0){
				this->inds0 = arrutilv2::copy_ptrarr_to_uptrarr(this->ncentroids, this->data_indices_init_from);
				tasks.push_back(this->set_C_Cl22s_from_inds0etc_ati());
				tasks.push_back(
					[this](TInt ti){
						if (ti == 0){

						}
					}
				);
			}
			
			else if (this->initialisation_method.compare("uniform") == 0){
				
				
				//tasks.emplace_back(	[this, ind0, ind1](TInt ti){ 
					//if (ti == 0){
						auto C_inds0 = this->get_C_inds0_uniform(ind0, ind1);
						this->C = std::move(std::get<0>(C_inds0));
						this->inds0 = std::move(std::get<1>(C_inds0));
						arrutilv2::set_rl22s(this->ncentroids, this->dimension, this->get_C(), this->get_C_l22s());
					//}
				//});
			}
			
			else if (this->initialisation_method.compare("kmeans++") == 0){
				
				this->inds0.reset(new TInt [this->ncentroids]);

				tasks.push_back( 
					[this, ind0, ind1](TInt ti){ 
						if (ti == 0){
							this->do_kmeanspp_initialisation(ind0, ind1);	
						}
					}
				);
			}
			
				
			
			else {
				throw std::runtime_error("expected initialisation method in {from_C, from_indices, uniform, kmeans++}. unrecognised initialisation scheme, bailing");
			}
			
			return tasks;
		}

		virtual TFloat get_validation_mse() = 0;
	
		//this is not exactly what we want for minibatch kmeans, but I will leave it here for now
		virtual void set_n_empty_clusters() final{
			this->n_empty_clusters = 0;
			for (TInt ci = 0; ci < this->ncentroids; ++ci){
				if (counts[ci] == 0){
					++this->n_empty_clusters;
				}
			}
		}


		virtual void set_X_tasks() = 0;
		/* the following are protected, as I don't want to hand pointers/(refs to non-const) out to the peanut gallery. Q: any advantage of this over just having these as protected variables ? A: (for unique_ptrs) not much. Derived classes can access raw pointer but not unique_ptr wrapper.*/

		TFloat * const get_valdata_l22s(){
			return valdata_l22s.get();
		}


		TFloat * const get_C(){
			return C.get();
		}				
		
		TFloat * const get_C_l22s(){
			return C_l22s.get();
		}

		TFloat * const get_sums(){
			return sums.get();
		}
		
		
		TInt * const get_counts(){
			return counts.get();
		}
		
		virtual void verbose_write_additional() = 0;
		
		
		
		/* function specific to minibatch versions,
		 *  use makeset_C_C_l22s_inds0_mati and set L of first batch (all others to 0)
		 * */
		 
		std::vector<std::function<void(TInt)> > minibatch_makeset_C_C_l22s_L_inds0_mati(minibatchapp::MiniBatchApp<TInt> & mba_in){
						
			std::vector<std::function<void(TInt)> > tasks;
			
			tasks = this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), mba_in.initialising_batch_size); //safe to move to densecentroid
			//set all L to "ncentroids" so that it is reflected that all change in zeroth round
			tasks.push_back(
				[this](TInt ti){
					if (ti == 0){
						std::fill_n(this->get_L(), this->ndata, this->ncentroids);
					}
				}
			);
			
			//set labels of initialising batch, update nchanges_on_batch[0]
			tasks.push_back(
				[this, &mba_in](TInt ti){
					TInt x0 = (ti*mba_in.initialising_batch_size)/this->nthreads;
					TInt x1 = ((ti+1)*mba_in.initialising_batch_size)/this->nthreads;
					
					this->set_L(x0, x1);
					//a hack:
					if (mba_in.nsubrounds > 1){
						mba_in.nchanges_on_batch[0] += x1 - x0;
					}
				}
			);
			
			return tasks;
		}
		
		
			//specific to minibatch kmeans (only some branches of this class are minibatch)
			std::function<void(TInt)> minibatch_subround_update(minibatchapp::MiniBatchApp<TInt> & mba_in){
			return [this, &mba_in] (TInt ti){
				if (ti == 0){
					mba_in.subround = (mba_in.subround + 1) % mba_in.nsubrounds;
					if (mba_in.subround == mba_in.nsubrounds - 1){
						this->nchanges = 0;
						for (auto & v : mba_in.nchanges_on_batch){
							this->nchanges += v;
							v = 0;
						}
						
					}
					else{
						this->nchanges = 10101;
					}
				}
			};
		}
			

		void minibatch_set_mse(const minibatchapp::MiniBatchApp<TInt> & mba_in){			
			if (this->round > 0 && (mba_in.subround == 0 || mba_in.subround == mba_in.nsubrounds - 1)){
				
				//this->round % mba_in.nsubrounds == 0){
				this->mse = this->getmeanl22at(); //note that this used labels determined with centroids other than the current ones. 
			}
			else{
				this->mse = -1;
			}
		}
		
		virtual TFloat getmeanl22at() = 0;
		
		
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_lowers_upper_mati(TInt data0, TInt data1){
			
			std::vector<std::function<void(TInt)> > tasks;
			tasks = this->base_makeset_C_C_l22s_inds0_mati(data0, data1);				
			tasks.push_back(
				[this, data0, data1](TInt ti){
					TInt ndata_loc = data1 - data0;
					TInt x0 = data0 + (ti*ndata_loc)/this->nthreads;
					TInt x1 = data0 + ((ti+1)*ndata_loc)/this->nthreads;
					this->set_upper_lowers_L(x0, x1);
					
				}
			);
			return tasks;
		}
		

		//only used by exact elkan kmeans. Currently for dense and sparse. A vickhale. 
		std::vector<std::function<void(TInt)> > exact_makeset_C_C_l22s_L_inds0_lowers_upper_mati(){
	
			std::vector<std::function<void(TInt)> > tasks = 
			this->makeset_C_C_l22s_L_inds0_lowers_upper_mati(static_cast<TInt> (0), this->ndata);
				
				
				//set mse. TODO : paralellise
				tasks.push_back(
				
				[this](TInt ti){	
					if (ti == 0 && this->get_initialisation_method().compare("kmeans++") != 0){
						this->mse = 0;
						for (TInt i = 0; i < this->getndata(); ++i){
							this->mse += (this->elkan_upper_base[i])*(this->elkan_upper_base[i]);
						}
						this->mse /= static_cast<TFloat>(this->ndata);
					}
					
				}
				
				);

			return tasks;
		}
		


		//will only be used by grow batch mse 
		std::function<void(TInt)> set_mse_sse_by_cluster_ati(const growbatchapp::GBApp<TInt, TFloat> & gbax, growbatchapp::GBMseApp<TInt, TFloat> & gbmseappx){
			
			
			return [this, &gbax, &gbmseappx](TInt ti){
				//minor computation, no point in parallelising
				if (ti == 0){
				
					std::fill(gbmseappx.sse_by_cluster.begin(), gbmseappx.sse_by_cluster.end(), 0);
					std::fill(gbmseappx.mse_by_cluster.begin(), gbmseappx.mse_by_cluster.end(), 0);
					
					for (TInt i = 0; i < gbax.ndata_active; ++i){
						gbmseappx.sse_by_cluster[this->L[i]] += gbmseappx.dn[i]*gbmseappx.dn[i];
					}
					
					
					for (TInt ci = 0; ci < this->ncentroids; ++ci){
						if (this->counts[ci] < 1){
							gbmseappx.mse_by_cluster[ci] = 0;
						}
						else{
							gbmseappx.mse_by_cluster[ci] = gbmseappx.sse_by_cluster[ci] / (this->counts[ci] - 1.);
						}
					}
				}
			};
		}



		

		//used by grow bath functions 
		std::vector<std::function<void(TInt)> > bgbmse_makeset_C_C_l22s_L_dn_inds0_mati(const growbatchapp::GBApp<TInt, TFloat> & gba){
			std::vector<std::function<void(TInt)> > tasks;
			tasks = this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), gba.ndata_active);
			/* set all L to "ncentroids" so that it is reflected that all change in zeroth round		*/
			tasks.push_back(
				[this](TInt ti){
					if (ti == 0){
						std::fill_n(this->get_L(), this->ndata, this->ncentroids);
						//std::fill_n(dn, this->ndata, 1e6);
					}
				}
			);
			//set labels and distances on active data
			tasks.push_back(
				this->set_L_dn_ati(0, gba.ndata_active)
			);
			return tasks;
		}
		
		
		//used by grow elkan bath functions 
		std::vector<std::function<void(TInt)> > bgbmse_makeset_C_C_l22s_L_lowers_dn_inds0_mati(const growbatchapp::GBApp<TInt, TFloat> & gba){
			std::vector<std::function<void(TInt)> > tasks;
			tasks = this->base_makeset_C_C_l22s_inds0_mati(static_cast<TInt> (0), gba.ndata_active);
			/* set all L to "ncentroids" so that it is reflected that all change in zeroth round		*/
			tasks.push_back(
				[this](TInt ti){
					if (ti == 0){
						std::fill_n(this->get_L(), this->ndata, this->ncentroids);
						//std::fill_n(dn, this->ndata, 1e6);
					}
				}
			);
			//set labels and distances on active data
			tasks.push_back(
				this->set_L_lowers_dn_ati(0, gba.ndata_active)
			);
			return tasks;
		}
		

		void gbmse_set_mse(const growbatchapp::GBApp<TInt, TFloat> & gba, const growbatchapp::GBMseApp<TInt, TFloat> & gbmseapp){
			if (gba.ndata_active != this->ndata){
				this->mse = -1;
			}
			else{
				this->mse = 0;
				for (TInt i = 0; i < this->ndata; ++i){
					this->mse += gbmseapp.dn[i]*gbmseapp.dn[i];
				}
				this->mse /= static_cast<TFloat> (this->ndata);
			}
		}
		
		

	public:

		virtual ~BaseDenseCentroidKmeans(){}

		BaseDenseCentroidKmeans(
		TInt nthreads, //1 
		TInt ndata, //2
		TInt dimension, //3
		TInt ncentroids, //4
		int cout_verbosity, //5
		int file_verbosity, //6
		std::ofstream & file, //7
		const std::string & initialisation_method, //8
		const TInt * const data_indices_init_from, //9
		bool setseed, //10
		TInt seed, //11
		TFloat maxtime, //12
		TInt maxrounds, //13 
		const std::string & verbose_filename, //14
		TInt nvaldata, //15
		TInt valperiod, //16
		const std::string & cmsewritefn, //17
		TInt cmserate, //18
		TFloat gbphi, //19
		const TFloat * const C_init //20 **** stays here, all others peal off and go to basecluster.

 
		):		
		
		cluster::BaseCluster<TInt, TFloat> (
		nthreads, //1
		ndata, //2
		dimension, //3
		ncentroids, //4 
		cout_verbosity, //5
		file_verbosity, //6
		file, //7
		initialisation_method, //8
		data_indices_init_from, //9
		setseed, //10
		seed, //11
		maxtime, //12
		maxrounds, //13
		verbose_filename, //14
		nvaldata, //15
		valperiod, //16
		cmsewritefn, //17
		cmserate, //18
		gbphi //19
		), //cmsewritefn, dmserate), 
		C_init(C_init),  data_l22s(new TFloat [ndata]), C(new TFloat [ncentroids*dimension]), C_l22s(new TFloat [ncentroids]), sums(new TFloat [ncentroids*dimension]), counts(new TInt [ncentroids]){
			
			
			valdata_l22s.reset(new TFloat [nvaldata] );
			
			
			this->checkinitialisationvalidity(this->C_init != nullptr);
			std::fill_n(this->get_sums(), this->ncentroids*this->dimension, 0);
			std::fill_n(this->get_counts(), this->ncentroids, 0);
			
		}
		
		/* do clustering and return : C, L, this->inds0 (this->inds0 may be nullptr, depending on initialisation method)
		 * duration, this->round, mse */
		std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> , TInt, TInt, TFloat> 
		get_6(){				
			
			//* initialise tasks */
			this->settasks();
			auto tstart = std::chrono::high_resolution_clock::now();
			/* run tasks */
						
			this->gopllcluster();
			
			auto tend = std::chrono::high_resolution_clock::now();
			TInt innerduration = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
			return std::make_tuple(std::move(C), std::move(this->L), std::move(this->inds0), innerduration, this->round, this->mse);
		}
		
		
		
		
		std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> > 
		/* C, LO, this->inds0 */
		get(){
			std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat> tup6 = this->get_6();
			return std::make_tuple(std::move(std::get<0>(tup6)), std::move(std::get<1>(tup6)), std::move(std::get<2>(tup6)));
		}
		
		std::tuple<TInt, TInt, TFloat >
		get_endstats(){
			auto tup6 = this->get_6();
			return std::make_tuple(std::get<3>(tup6), std::get<4>(tup6), std::get<5>(tup6));
		}
		
		
		
		const TFloat * const get_data_l22s(){
			return data_l22s.get();
		}

		const TFloat * const get_C_init(){
			return C_init;
		}

		//Return an estimate of memory requirement of this base class
		virtual TInt get_approximate_memory_requirement() = 0;
		
		
		void EB_verbose_write_additional(){
		
			this->get_verbose_file() << "\nlowers base:\n";
			for (TInt ci = 0; ci < this->getncentroids(); ++ci){
				this->get_verbose_file() << this->elkan_lowers_base[ci] << "\t";
			}
			this->get_verbose_file() << "\n\nupper base:\n" << this->elkan_upper_base[0] << "\n";
		}
		
};

}

//extern template class kmeans::BaseDenseCentroidKmeans<size_t, double>;
//extern template class kmeans::BaseDenseCentroidKmeans<size_t, float>;





#endif



		//std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> , TInt, TInt, TFloat>
		//get_6_memexcept(){
			
			//arrutilv2::proxy_openblas_set_num_threads(1);
			//std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat> tup;
			
			
			//std::string ENTERING_OPENBLAS_NUM_THREADS = "";
			//if (std::getenv("OPENBLAS_NUM_THREADS")){
				//ENTERING_OPENBLAS_NUM_THREADS = std::getenv("OPENBLAS_NUM_THREADS");
			//}

			//try{
				//retyrb get_6();
			//}
			
			//catch (std::bad_alloc& ba){
			
				//setenv("OPENBLAS_NUM_THREADS", ENTERING_OPENBLAS_NUM_THREADS.c_str(),1);
				//std::cerr << "bad_alloc caught: " << ba.what() << '\n';
				//std::cerr << "the attempt to cluster failed. Probably due to a memory allocation failure : the memory demand was too large (ndata/ncentroids/ngroups/n... too large). For now, will return tuple of nullptrs etc\n" << std::endl; 
			//}
		
		//return 	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat > {nullptr, nullptr, nullptr, 0, 0, 0};
	//}
			
		//}

