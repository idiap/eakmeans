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

#ifndef PLL_BASECLUSTER_H
#define PLL_BASECLUSTER_H

#include <memory>
#include <iostream>
#include <vector>
#include <mutex>
#include <atomic>
#include <thread>
#include <tuple>
#include <limits>

#include "barrierutil.h"
#include "arrutilv2l3.h"
#include "arrutilv2mse.h"
#include "arrutilv2discrete.h"

#include "sample.h"
#include "initialise2.h"
#include "stringutilclustering.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

namespace cluster{

/* Considered making it like  : BInt, SInt, TFloat -- size_t, unsigned short, double.
 * and using SInt for L, ncentroids. But it is not clear that it will speed things up, due to the "placenta"
 * Before undertaking this potentially large change it would be worth performing experiments determing speeds of multiplying, adding various int types with other int types.  */
template <typename TInt, typename TFloat>
class BaseCluster{
	
	private:
		
		void checknthreadsverbositycompatiblity(TInt nthreads, int file_verbosity, const std::string & verbose_filename){
			if (file_verbosity > 1 && nthreads > 1){
				std::string strangerequest("performance warning : extreme verbosity to file (file_vebosity > 1) is enabled for multithreading, but note that extreme file verbosity is only intended to be used for creating animations, and so it may slow down performace with multithreading. That said, it should work (file verbosity processed in thread 0).");
				std::cout << strangerequest << std::endl;
			}
			
			if (file_verbosity <= 1 && verbose_filename.compare("") != 0){
				std::string res = std::string("Received non-extreme file_verbosity value : ")  + std::to_string(file_verbosity) + " with non-empty verbose_filename : " + verbose_filename + ". This is confusing as verbose_filename is only relevant when file_verbosity >= 2. To prevent confusion, it is required that if file_verbosity <= 1, then verbose_filename is the empty string.";
				throw std::runtime_error(res);
			}
		}
		
		
		void checkvalidationcompatibility(int cver, int fver, TInt nthreads, bool withvalidation){
			if (withvalidation) {
				if (!(cver == 2 || fver > 0)){
					std::string err = std::string("Request to use validation set, with incompatible verbosity settings. To run with a validation set, eith cver == 2 or fvaer >0. Currently, cver = ") + std::to_string(cver) + " and fver = " + std::to_string(fver) + ".";
					throw std::runtime_error(err);
				}
				if (nthreads > 1){
					std::cerr << "Request to run with a validation set using multiple (" << std::to_string(nthreads) << ") threads. Request will be processed, although validation set will be processed on a single threads. The use of a validation set may considerably slow down k-means, as there is no optimisation in the assignment for the validation set. Use of the validation set should only be used for illustatitive/comparative purposes.  ";
				}
			}
		}
		
		
		std::vector<std::function<void(TInt) > > section_tasks;
		
		virtual void set_section_tasks(){
						
			set_C_tasks();
			section_tasks = C_tasks;

			set_X_tasks();
			section_tasks.insert(section_tasks.end(), X_tasks.begin(), X_tasks.end());
		
		}
				
		void updateduration(){
			auto tcurrent = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::milliseconds>(tcurrent - tstart).count();
		}
	
	
		void endroundupdate(){
			this->iscomplete = (nchanges == 0) || (duration > maxtime) || (round >= maxrounds);
			this->nchanges = 0;
			++this->round;
		}


		virtual void set_summaries() = 0;
		
		//TODO : change mse to energy, as may be energy other than l2
		virtual void set_mse() = 0;



	
	
	protected:

		//Certain of these variables `should' be private with getters 
		std::string algname;
		TInt nthreads;
		TInt ndata;
		TInt dimension;
		TInt ncentroids; 
		int cout_verbosity;
		int file_verbosity; 
		std::ofstream * fileptr;
		std::string initialisation_method;
		const TInt * data_indices_init_from;
		bool setseed;
		TInt seed;
		TFloat maxtime;
		TInt maxrounds;
		std::ofstream verbose_file;
		TInt valperiod;
		TInt nvaldata;
		TInt nchanges;
		std::minstd_rand0 g;
		std::unique_ptr<TInt []> L;
		bool withvalidation;
		std::atomic<TInt> ndcalcs_notX;
		std::atomic<TInt> ndcalcs_X;
		TInt round;
		std::mutex work_mutex;
		std::unique_ptr<TInt []> inds0;
		TFloat mse;
		std::chrono::time_point<std::chrono::high_resolution_clock> tstart;
		TFloat duration;
		TInt approximate_memory_requirement;
		bool iscomplete;
		TInt n_empty_clusters;
		TFloat val_mse;

		//printing actions to std out
		std::function<void()> cout_startsummary;
		std::function<void()> cout_roundsummary;
		std::function<void()> cout_finalsummary;

		//printing actions to file
		std::function<void()> file_startsummary;
		std::function<void()> file_roundsummary;
		std::function<void()> file_finalsummary;

		//manage printing actions to stdout and to file
		std::function<void()> startsummary;
		std::function<void()> roundsummary;
		std::function<void()> finalsummary;

	
	virtual TFloat get_validation_mse() = 0;	


	/* {cout_startsummary, file_startsummary} ---> 
	 * startsummary, etc controlling when mse and n_empty_clusters are set.
	 * */
	void combine_cout_file_summaries(){
		
		startsummary = [this](){
			if (withvalidation &&  (cout_verbosity > 0 || file_verbosity > 0)){
				this->val_mse = this->get_validation_mse();
			}
			else{
				this->val_mse = -1;
			}
			cout_startsummary();
			file_startsummary();
		};
		roundsummary = [this](){
			
			if (cout_verbosity > 1 || file_verbosity > 0){
				this->set_mse();
				if (withvalidation && this->round % valperiod == 0){
					this->val_mse = this->get_validation_mse();
				}
				else{
					this->val_mse = -1;
				}
			}
			
			cout_roundsummary(); 
			file_roundsummary();
		};
		finalsummary = [this](){
			//only if verbosity is non-(00) will mse and n_empty_clusters be used
			if (cout_verbosity > 0 || file_verbosity > 0){
				this->set_n_empty_clusters();
				this->set_mse();
			}
			
			if (withvalidation & (cout_verbosity > 0 || file_verbosity > 0)){
				this->val_mse = this->get_validation_mse();
			}
			
			else{
				this->val_mse = -1;
			}
			cout_finalsummary();
			file_finalsummary();
		};
	}

	

	
	
	//input parameter should be true if the variable(s) necessary for a C automatic initialisation are provided to constructor of derived
	void checkinitialisationvalidity(bool with_C_initialisation){
		
		if (with_C_initialisation == false){					
			if (this->data_indices_init_from != nullptr && this->initialisation_method.compare("from_indices") != 0){
				throw std::logic_error("non-nullptr this->data_indices_init_from and yet this->initialisation_method != from_indices, bailing");
			}
			else if (this->data_indices_init_from == nullptr && this->initialisation_method.compare("from_indices") == 0){
				throw std::logic_error("this->data_indices_init_from is nullptr and yet this->initialisation_method == from_indices, bailing");
			}
		}
		
		else{
			if (this->data_indices_init_from != nullptr){
				throw std::logic_error("this->data_indices_init_from is non-nullptr and C_init is non-nullptr, bailing");
			}
			else if (this->initialisation_method.compare("from_C") != 0){
				throw std::logic_error("C_init is non-nullptr and this->initialisation_method != from_C, bailing");
			}
		}
	}



	void gopllcluster(){

		/*set approximate_memory_requirement, not done in constructor as would not call derived virtual function*/
		approximate_memory_requirement = get_approximate_memory_requirement();	
		
		//set summaries (ditto above comment)
		set_summaries();

		tstart = std::chrono::high_resolution_clock::now();
		duration = 0;

		auto ninittasks = initialisation_tasks.size();		
		auto nsections = section_tasks.size();

		std::vector<std::function<void()>> initialisationend_tasks (ninittasks, [](){});
		std::vector<std::function<void()>> sectionend_tasks (nsections, [](){}); 
		std::function<void()> closing_task = finalsummary;//[](){};

	


		initialisationend_tasks[ninittasks -1] = [this](){  
			startsummary();	
			updateduration();
			roundsummary();
			nchanges = 0;
			++round;
		};

			 
		sectionend_tasks[nsections - 1] = [this](){
			this->updateduration();
			this->roundsummary();
			this->endroundupdate();
		};
		
		//closing_task = finalsummary;
			
		/* perform clustering */
		
		
		

		stdthreadutil::launch_btasks_rbtasks(
		nthreads, 
		initialisation_tasks, initialisationend_tasks, 
		section_tasks, sectionend_tasks, [this](){return this->iscomplete;},
		closing_task);
	
	}
			
		virtual void verbose_write_additional(){};
		std::vector<std::function<void(TInt) > > initialisation_tasks;
		std::vector<std::function<void(TInt) > > C_tasks;
		std::vector<std::function<void(TInt) > > X_tasks;

		virtual void set_initialisation_tasks() = 0;

		//the decomposition into C tasks and X tasks is arbitrary, it is used only to help me understand my code. It is possible to put an X task in C tasks, and vica versa.
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;

		void settasks(){
			/* set the serious tasks */
			set_initialisation_tasks();
						
			set_section_tasks();

			/* if extreme verbosity, add additional writing task (save centroids etc in each round) */
			if (file_verbosity > 1){
				
				initialisation_tasks.emplace_back(
					[this](TInt ti){
						if (ti == 0){
							verbose_write();
						}
					} 
				);
				
				section_tasks.emplace_back(
					[this](TInt ti){
						if (ti == 0){
							verbose_write();
						}
					} 
				);
			}
		}
		
		void setalgname(std::string newalgname){
			this->algname = newalgname;
		}

				
		std::ofstream & get_verbose_file(){
			return verbose_file;
		}
		
		void set_silent_summaries(){
			//default for 0 verbosity			
			cout_startsummary = [](){};
			cout_roundsummary = [](){};
			cout_finalsummary = [](){};
	
			//default for 0 verbosity
			file_startsummary = [](){};
			file_roundsummary = [](){};
			file_finalsummary = [](){};
		}	
	
	
		void set_summaries_exact(){//certain derived classes will use this function, but as they live on different branches I place this function here. I know that this goes against the tried and tested "is-a" approach, buy can't see a cleaner solution.
		
		
			set_silent_summaries();
	
			if (cout_verbosity == 1){
				cout_startsummary = [this](){
					std::cout << 
					stringutil::clustering::pll::exact::getstartsummary_v1(algname, approximate_memory_requirement, mse, val_mse)
					 << std::flush;
				};
				
				cout_roundsummary = [this](){
					std::cout << 
					stringutil::clustering::pll::exact::getroundsummary_v1(nchanges) 
					<< std::flush;
				};
	
				cout_finalsummary = [this](){
					std::cout << 
					stringutil::clustering::pll::exact::getfinalsummary_v1(round - 1, nchanges, static_cast<TInt>(ndcalcs_X), static_cast<TInt>(ndcalcs_notX + ndcalcs_X), duration, mse, n_empty_clusters, val_mse) 
					<< std::flush;
				};
			}
		
			else if (cout_verbosity == 2){
				cout_startsummary = [this](){
					std::cout <<
					stringutil::clustering::pll::exact::getstartsummary_v2(algname, approximate_memory_requirement, mse, val_mse)
					<< std::flush;
				};
				
				cout_roundsummary = [this](){
					std::cout << 
					stringutil::clustering::pll::exact::getroundsummary_v2(round, nchanges, static_cast<TInt>(ndcalcs_X), static_cast<TInt>(ndcalcs_notX + ndcalcs_X), duration, mse, val_mse)	
					 << std::flush;
				};
	
				cout_finalsummary = [this](){
					std::cout << stringutil::clustering::pll::exact::getfinalsummary_v2(round - 1, nchanges, static_cast<TInt>(ndcalcs_X), static_cast<TInt>(ndcalcs_notX + ndcalcs_X), duration, mse, n_empty_clusters, val_mse)
					<< std::flush;
				};
			}
						
			
			
			if (file_verbosity == 1){
				 throw std::runtime_error("\nI need to implement file_verbosity == 1 to have sensible feedback here \n ");
			}
			
			if (file_verbosity == 2 || file_verbosity == 3){
				file_startsummary = [this](){
					*fileptr << 
					stringutil::clustering::pll::exact::getstartsummary_v2(algname, approximate_memory_requirement, mse, val_mse)
					 << std::flush;
				};
				
				file_roundsummary = [this](){
					*fileptr << 
					stringutil::clustering::pll::exact::getroundsummary_v2(round, nchanges,  static_cast<TInt>(ndcalcs_X), static_cast<TInt>(ndcalcs_notX + ndcalcs_X), duration, mse, val_mse)
					 << std::flush;
				};
	
				file_finalsummary = [this](){
					*fileptr << 
					stringutil::clustering::pll::exact::getfinalsummary_v2(round - 1, nchanges, static_cast<TInt>(ndcalcs_X), static_cast<TInt>(ndcalcs_notX + ndcalcs_X), duration, mse, n_empty_clusters, val_mse)
					 << std::flush;
				};
			}
			
			combine_cout_file_summaries();
			
		}
			
						

	
	public:

	BaseCluster(TInt nthreads, TInt ndata, TInt dimension, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const std::string & initialisation_method, const TInt * const data_indices_init_from, bool setseed, TInt seed, TFloat maxtime, TInt maxrounds, const std::string & verbose_filename, TInt valperiod, TInt nvaldata):
		algname("virtual cluster base"), nthreads(nthreads), ndata(ndata), dimension(dimension), ncentroids(ncentroids), cout_verbosity(cout_verbosity), file_verbosity(file_verbosity), fileptr(&file), initialisation_method(initialisation_method), data_indices_init_from(data_indices_init_from), setseed(setseed), seed(seed), maxtime(maxtime), maxrounds(maxrounds), valperiod(valperiod), nvaldata(nvaldata), nchanges(ndata), L(new TInt [ndata]), ndcalcs_notX(0), ndcalcs_X(0), round(0), iscomplete(false)
		{
			withvalidation = (nvaldata != 0);
			checknthreadsverbositycompatiblity(nthreads, file_verbosity, verbose_filename);
			checkvalidationcompatibility(cout_verbosity, file_verbosity, nthreads, withvalidation);

			if (setseed == true){
				this->seed = seed;
			}
			
			else{
				this->seed = time(NULL);
			}
			
			g.seed(this->seed);
			srand(this->seed);
			

			if (file_verbosity > 1){
				verbose_file.open(verbose_filename, std::ofstream::app);
				//TODO : there is no file clearing before writing, (in the case of the hulk, done externally)
				verbose_file << "\n|||||||||||||||||||||\nverbose output initialisation from  " << algname << "\n\n\n";
			}
		
			startsummary = [](){};
			roundsummary = [](){};
			finalsummary = [](){};

		}
		
		

		void verbose_write(){
			verbose_write_additional();
		}
				


		const std::string & getalgname(){
			return algname;
		}
		
		
		TInt getnthreads(){
			return nthreads;
		}
		
		TInt getndata(){
			return ndata;
		}
		
		TInt getdimension(){
			return dimension;
		} 
		
		
		TInt getncentroids(){
			return ncentroids;
		}
		
		int get_cout_verbosity(){
			return cout_verbosity;
		}
		
		int get_file_verbosity(){
			return file_verbosity;
		}
		
		std::ofstream * const getfileptr(){
			return fileptr;
		}
		
		const std::string & get_initialisation_method(){
			return initialisation_method;
		}
		
		
		const TInt * const get_data_indices_init_from(){
			return data_indices_init_from;
		}
		
		
		TInt getseed(){
			return seed;
		}
		
		TFloat getmaxtime(){
			return maxtime;
		}
		
		TInt getmaxrounds(){
			return maxrounds;
		}
		
		TInt * const get_L(){
			return L.get();
		}

		
		virtual TInt get_approximate_memory_requirement() = 0;


		virtual void set_n_empty_clusters() = 0;
		



		virtual ~BaseCluster(){
			if (verbose_file.is_open()){
				verbose_file.flush();
				verbose_file.close();
			}
		}
	};
} 
 


//extern template class cluster::BaseCluster<size_t, double>;
//extern template class cluster::BaseCluster<size_t, float>;



#endif
