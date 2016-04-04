#ifndef PLL_BASECLUSTER_H
#define PLL_BASECLUSTER_H

#include <memory>
#include <iostream>
#include <vector>
#include <mutex>
#include <atomic>
#include <thread>
#include <chrono>
#include <tuple>
#include <limits>

#include "barrierutil.h"
#include "stringutilclustering.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>


#include "minibatchapp.h"
#include "growbatchapp.h"


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

			section_tasks = {};


		


						
			set_C_tasks();
			section_tasks.insert(section_tasks.end(), C_tasks.begin(), C_tasks.end());
			//section_tasks = C_tasks;
			

			set_X_tasks();
			section_tasks.insert(section_tasks.end(), X_tasks.begin(), X_tasks.end());
			
			
			
			if (this->writecleanmses == true){
				section_tasks.push_back(
					[this](TInt ti){
						if (ti == 0){
							auto time_now = std::chrono::high_resolution_clock::now();
							TInt duration_since_last = std::chrono::duration_cast<std::chrono::milliseconds>(time_now - t_last_stats_gather).count();
							
							//if (less than a 2 seconds have passed since start / long time since last / no changes
							if ((duration < 2000) || (duration_since_last > this->cmserate) || this->nchanges == 0){
							//only do if times so desire, or at end (no changes)
								updateduration();								
								auto clean_mse = this->base_get_clean_mse();								
								this->cmsefile << duration << " \t" << clean_mse << "\n";
								this->cmsefile.flush();								  								
								t_last_stats_gather = std::chrono::high_resolution_clock::now();
							}
							

						}
					}
				);	
			}

		}
				
		void updateduration(){
			auto tcurrent = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::milliseconds>(tcurrent - tstart).count();
			duration -= gathering_statistics_duration;
		}
	
	
		virtual void endroundupdate(){ //this works for most algorithms, but not for growbatches.
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
		TInt nvaldata;
		TInt valperiod;
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
		//Time doing meaningless statistics gathering, should be subtracted.
		//Note that this hack should only be used in base_get_clean_mse
		TFloat gathering_statistics_duration; 
		//Total time excluding gathering_statistics_duration. 
		std::chrono::time_point<std::chrono::high_resolution_clock> t_last_stats_gather;
		TFloat duration; 
		TInt approximate_memory_requirement;
		bool iscomplete;
		TInt n_empty_clusters;
		TFloat val_mse;
		std::string cmsewritefn;
		std::ofstream cmsefile;
		TInt cmserate;
		TFloat gbphi;
		bool writecleanmses;

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

	
	
	//compute mse from scratch

	virtual TFloat get_clean_mse(bool on_validation){
		throw std::logic_error("get clean mse called in basecluster makes no sense, should be call to overriden version (sparse/dense)");
		//using namespace std::literals;
		//std::this_thread::sleep_for(std::chrono::seconds(1));
		
	};


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
	
	//Stop the clock and get a clean mse
	TFloat base_get_clean_mse(){
		std::chrono::time_point<std::chrono::high_resolution_clock> tstart = std::chrono::high_resolution_clock::now();
		bool cmse_on_validation = (this->nvaldata != 0);
		TFloat clean_mse = this->get_clean_mse(cmse_on_validation);
		std::chrono::time_point<std::chrono::high_resolution_clock> tend = std::chrono::high_resolution_clock::now();
		this->gathering_statistics_duration += 
		std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
		return clean_mse;
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
		gathering_statistics_duration = 0;
		
		auto ninittasks = all_initialisation_tasks.size();		
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
		all_initialisation_tasks, initialisationend_tasks, 
		section_tasks, sectionend_tasks, [this](){return this->iscomplete;},
		closing_task);
	
	}
			
		virtual void verbose_write_additional(){};
		std::vector<std::function<void(TInt) > > all_initialisation_tasks; //will include l22s in basedensecentroidkmeans
		std::vector<std::function<void(TInt) > > initialisation_tasks;
		std::vector<std::function<void(TInt) > > C_tasks;
		std::vector<std::function<void(TInt) > > X_tasks;

		virtual void set_all_initialisation_tasks() = 0;

		//the decomposition into C tasks and X tasks is arbitrary, it is used only to help me understand my code. It is possible to put an X task in C tasks, and vica versa.
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;

		void settasks(){
			/* set the serious tasks */
			set_all_initialisation_tasks();
						
			set_section_tasks();

			/* if extreme verbosity, add additional writing task (save centroids etc in each round) */
			if (file_verbosity > 1){
				
				all_initialisation_tasks.emplace_back(
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
	
		//certain derived classes will use this function, but as they live on different branches I place this function here. I know that this goes against the tried and tested "is-a" approach, buy can't see a cleaner solution.
		void set_summaries_exact(){
		
		
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
		
	
	
	//as with set_summaries_exact, only certain derived classes will use this function, but as they live on different branches I place this function here. Considered solution with friend function, in the end with this. There is probably a neater solution, but would require major refactoring. 	
	void set_summaries_minibatch(const minibatchapp::MiniBatchApp<TInt> & mba){

		
		this->set_silent_summaries();
		
		if (this->cout_verbosity == 1){
			this->cout_startsummary = [this](){   
				std::cout << 
				stringutil::clustering::pll::minibatch::getstartsummary_v1(this->algname, this->approximate_memory_requirement, this->val_mse)
				 << std::flush;
			};
			
			this->cout_roundsummary = [this, &mba](){   
				std::cout << 
				stringutil::clustering::pll::minibatch::getroundsummary_v1(mba.subround, mba.nsubrounds, this->nchanges) 
				<< std::flush;
			};
	
			this->cout_finalsummary = [this, &mba](){   
				std::cout << 
				stringutil::clustering::pll::minibatch::getfinalsummary_v1(this->round, mba.nsubrounds, this->nchanges, static_cast<TInt>(this->ndcalcs_X), static_cast<TInt>(this->ndcalcs_notX + this->ndcalcs_X), this->duration, this->mse, this->n_empty_clusters, this->val_mse) 
				<< std::flush;
			};
		}
	
		else if (this->cout_verbosity == 2){
			
			this->cout_startsummary = [this](){   
				std::cout <<
				stringutil::clustering::pll::minibatch::getstartsummary_v2(this->algname, this->approximate_memory_requirement, this->val_mse)
				<< std::flush;
			};
			
			this->cout_roundsummary = [this, &mba](){  
				std::cout << 
				stringutil::clustering::pll::minibatch::getroundsummary_v2(this->round/mba.nsubrounds, mba.nsubrounds, mba.subround, this->nchanges, static_cast<TInt>(this->ndcalcs_X), static_cast<TInt>(this->ndcalcs_notX + this->ndcalcs_X), this->duration, this->mse, this->val_mse)	
				 << std::flush;
				 
			};
	
			this->cout_finalsummary = [this, &mba](){   
				std::cout << stringutil::clustering::pll::minibatch::getfinalsummary_v2(this->round, mba.nsubrounds, this->nchanges, static_cast<TInt>(this->ndcalcs_X), static_cast<TInt>(this->ndcalcs_notX + this->ndcalcs_X), this->duration, this->mse, this->n_empty_clusters, this->val_mse)
				<< std::flush;
			};
		}
					
		this->combine_cout_file_summaries();
	}
	
	//as with set_summaries_minibatch, we include this function here:

	void set_summaries_growbatch(const growbatchapp::GBApp<TInt, TFloat> & gbapp)
		{
			
			this->set_silent_summaries();
			
			if (this->cout_verbosity == 1){
				this->cout_startsummary = [this](){ 
					std::cout << 
					stringutil::clustering::pll::growbatch::getstartsummary_v1(this->algname, this->approximate_memory_requirement, this->val_mse)
					 << std::flush;
				};
				
				this->cout_roundsummary = [this, &gbapp](){ 
					std::cout << 
					stringutil::clustering::pll::growbatch::getroundsummary_v1(this->nchanges, gbapp.ndata_active != gbapp.ndata_active_previous) 
					<< std::flush;
				};
	
				this->cout_finalsummary = [this](){ 
					std::cout << 
					stringutil::clustering::pll::growbatch::getfinalsummary_v1(this->round, static_cast<TInt>(this->ndcalcs_X), static_cast<TInt>(this->ndcalcs_notX + this->ndcalcs_X), this->duration, this->mse, this->n_empty_clusters, this->val_mse) 
					<< std::flush;
				};
			}
		
			else if (this->cout_verbosity == 2){

				this->cout_startsummary = [this](){ 
					std::cout <<
					stringutil::clustering::pll::growbatch::getstartsummary_v2(this->algname, this->approximate_memory_requirement, this->val_mse)
					<< std::flush;
				};
				
				this->cout_roundsummary = [this, &gbapp](){
					std::cout << 
					
			
					
					stringutil::clustering::pll::growbatch::getroundsummary_v2(this->round, gbapp.ndata_active, gbapp.d_C__over__d_AB, this->nchanges, static_cast<TInt>(this->ndcalcs_X), static_cast<TInt>(this->ndcalcs_notX + this->ndcalcs_X), this->duration, this->mse, this->val_mse)	
					 << std::flush;
					 
				};
	
	

				this->cout_finalsummary = [this](){ 
					std::cout << stringutil::clustering::pll::growbatch::getfinalsummary_v2(this->round, static_cast<TInt>(this->ndcalcs_X), static_cast<TInt>(this->ndcalcs_notX + this->ndcalcs_X), this->duration, this->mse, this->n_empty_clusters, this->val_mse)
					<< std::flush;
				};
			}
						
			this->combine_cout_file_summaries();

		}			
			

		//used by growbatch classes (dense and sparse)
		void BGB_constructor_helper(const TInt & batchsize0, growbatchapp::GBApp<TInt, TFloat> & gbapp){ //gbapp should be be passed as this->gba, it will be set-up here

			//TODO : I have changed this to batchsize0 from batchsize0, which might break Partion Split Kmeans. Sort this out when have time, now writing a paper...
			gbapp.ndata_active = batchsize0; //std::max(this->ndata/100, 2*batchsize0);

			if (gbapp.ndata_active > this->ndata){
				std::string rer("The initial batch size of ");
				rer = rer + std::to_string(batchsize0) + 
				std::string(" is too great. Consider a smaller initial batch size, one which is smaller than half the amount of data (") + 
				std::to_string(this->ndata) + 
				std::string(").");
				
				throw std::runtime_error(rer);
			}
				
			
			gbapp.ndata_active_previous = 0;
			
			this->setalgname("Base Grow Batch");
			gbapp.growthfactor = 2.0;
			gbapp.threshold = 1.0; //TODO : is threshold humpty dumpty?
			gbapp.d_C__over__d_AB = 42.0;		
					
			
			gbapp.delta_C.reset(new TFloat [this->ncentroids]);
			
			
			if (this->ndata % 2 != 0){
				this->ndata -= 1;
			}
		}		

		

		//should only be used by growbatch kmeans
		virtual bool should_double(){
			throw std::logic_error("should_double called at base (basecluster) : either the class calling it should not, or it should override this the base implementation");
			return false;
		}
			

		//will only be used by growbatch kmeans
		std::function<void(TInt ti)> update_ndata_active_ati(growbatchapp::GBApp<TInt, TFloat> & gbapp){
		
		return 	
		[this, &gbapp](TInt ti){
			if (ti == 0){
				
				TInt outgoing_ndata_active = gbapp.ndata_active;
				
				//std::cout << "ndata_active : " << this->gba.ndata_active  << "\t delta_C : " << delta_C2_scalar << "\t delta_C2_AB : " << delta_C2_AB_scalar << std::endl;
				
				
				
				if (this->should_double() == true){
					if (gbapp.ndata_active == this->ndata){
						//std::cout << "Should stop : ndata_active == ndata and delta_C2_AB_scalar > threshold*delta_C2_scalar" << std::endl;
					}
					
					else{
						if (gbapp.ndata_active * 2 >= this->ndata){
							//std::cout << "`Growing active to full data'. " << std::endl;
							gbapp.ndata_active = this->ndata;
						}
						
						else{
							//std::cout << "`end of C : Doubling'. " << std::endl;
							gbapp.ndata_active *= 2;
						}
					}
				}
				
				gbapp.ndata_active_previous = outgoing_ndata_active;

			}
		};
	}


	//used only by grow batch algorithms
	bool gbmse_should_double(growbatchapp::GBApp<TInt, TFloat> & gbapp, const TFloat * const mse_by_cluster){

		
		/* ratios will be ((sqrt)mse_by_cluster / delta_C) for each centroid */
		std::vector<TFloat> ratios (this->ncentroids, 0);
		for (TInt ci = 0; ci < this->ncentroids; ++ci){
			if (gbapp.delta_C[ci] < 1e-10){
				/* The centroid did not move: large ratio for (mse / delta_C) */
				ratios[ci] = 1e13;
			}
			
			else if (mse_by_cluster[ci] < 1e-10){
				/* The mse is zero, this means that all data the same for this centroid. In the theory it could also be that there is no data for the centroid, but this case should be caught above. */
				ratios[ci] = 1e-13;
			}
			
			else{
				/* Note that large ratio means ==> rmse > delta C ==> more data should be added. */
				ratios[ci] = sqrt(mse_by_cluster[ci]) / gbapp.delta_C[ci]; 
			}
		}
		

		/* propotion of ratios required to vote for no increase to get no increase. [increase rank --> more doubling] */
		TFloat rank = 0.5; 
		/* delta_C > phi * rmse --> vote for no increase. [increase phi --> more doubling] */
		//gbphi : this is maybe a bit low. could experiment with at some point?
		std::nth_element(ratios.begin(), ratios.begin() + ratios.size()*rank, ratios.end());

		gbapp.d_C__over__d_AB = 1./ratios[ratios.size()*rank];
		
		
		/* Note d_C__over__d_AB < phi ==> 
		 * 1./ratios[ratios.size()*rank] < phi ==>
		 * ratios[ratios.size()*rank] >= 1/phi ==>
		 * less than (100 * rank)% of ratios less than 1/phi ==>
		 * less than (100 * rank)% of (rmse / delta_C) less than 1/phi ==> 
		 * more than (100 * (1 - rank))% of rmse / delta_C greater than 1/phi ==>
		 * more than (100 * (1 - rank))% of delta_C / rmse less than phi ==>
		 * (having delta_C /rmse < phi ==> delta_C < phi * rmse ==> vote for more data)
		 * more than (100 * (1 - rank))% of ratios vote for more data.
		 * less than (100 * rank)% vote for the same amount of data ==>
		 * GET MORE DATA!
		 * */


		//if (gbapp.ndata_active != gbapp.ndata_active_previous){
			//return false;
		//}
		



		return gbapp.d_C__over__d_AB < this->gbphi;
		

	}
	
	
	
		//will only be used by grow batch mse 
		void BGBM_constructor_helper(growbatchapp::GBMseApp<TInt, TFloat> & gbmseapp){
			gbmseapp.sse_by_cluster = std::vector<TFloat> (this->ncentroids, 0);			
			gbmseapp.mse_by_cluster = std::vector<TFloat> (this->ncentroids, 0);			
			gbmseapp.dn.reset(new TFloat [this->ndata]);
			this->setalgname("Dense/Sparse Base Grow Batch Mse (should be set in specific class)");
		}
		

						

	
	public:

	BaseCluster(
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
	TFloat gbphi = 0
	):
		algname("virtual cluster base"), nthreads(nthreads), ndata(ndata), dimension(dimension), ncentroids(ncentroids), cout_verbosity(cout_verbosity), file_verbosity(file_verbosity), fileptr(&file), initialisation_method(initialisation_method), data_indices_init_from(data_indices_init_from), setseed(setseed), seed(seed), maxtime(maxtime), maxrounds(maxrounds), nvaldata(nvaldata), valperiod(valperiod), nchanges(ndata), L(new TInt [ndata]), ndcalcs_notX(0), ndcalcs_X(0), round(0), iscomplete(false)
		{
			
			
			
			
			
			
			this->cmserate = cmserate;

			this->gbphi = gbphi;

			std::cout << "The value of gbphi is " << this->gbphi << std::endl;
			
			this->cmsewritefn = cmsewritefn;


			this->writecleanmses = false;
			if (this->cmsewritefn.compare("") != 0){
				
				this->writecleanmses = true;
				this->cmsefile.open(this->cmsewritefn);//, std::ofstream::ios);
				
			}

			withvalidation = (nvaldata != 0 &&  this->writecleanmses == false);
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
			
			this->cmsefile.close();
			if (verbose_file.is_open()){
				verbose_file.flush();
				verbose_file.close();
			}
		}
	};
} 
 

#endif


//extern template class cluster::BaseCluster<size_t, double>;
//extern template class cluster::BaseCluster<size_t, float>;
		////std::cout << "In gbmse_should_double... " << std::flush;
		////for (TInt q = 0; q < 100; ++q){
			////std::cout << rand()%100 << " " << std::flush << std::endl;
		////}
		//if (rand()%100 > 5){
			////std::cout << "Don't double! " << std::flush; 
			//return false;
		//}
		
		//else{
			
		//}
		//if (this->round %3 != 0){
			//return false;
		//}







				
				//std::cout << "cmsewritefn in basecluster set to : " << this->cmsewritefn <<  std::endl;
				//std::abort();

				//TODO: a check that the path is fiable. 
								//if (this->nchanges == 0){
									//std::cout << "Final cmsefile write with nchanges == 0 taking place" << std::endl;
								//}
