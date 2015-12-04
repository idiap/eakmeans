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

#ifndef PLL_BASEKMEANSSS_H
#define PLL_BASEKMEANSSS_H

#include "basecluster.h"

#include "arrutilv2mse.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseKmeans : public cluster::BaseCluster<TInt, TFloat> {

	//As this is a template class, inherited variables either need to be referenced as BaseCluster<TInt, TFloat>::v or as this->v or by declaring them up front as we do below. I have had issues with the later, and am so going with the -> option.
	
	
	private:
	
		virtual void set_mse() = 0;
		


		virtual void set_summaries() = 0;
		
			


	protected:
	
		//variables specific to dense l2:
		const TFloat * data;
		const TFloat * C_init;
		std::unique_ptr<TFloat []> data_l22s;
		std::unique_ptr<TFloat []> C;
		std::unique_ptr<TFloat []> C_l22s;
		std::unique_ptr<TFloat []> sums;
		std::unique_ptr<TInt []> counts; 
		std::unique_ptr<TFloat []> dsums;
		std::unique_ptr<TInt []> dcounts;
		const TFloat * valdata;
		std::unique_ptr<TFloat []> valdata_l22s;
		TInt maxpermultiplyblock_for_validation;


		std::function<void(TInt)> base_set_S_H_ati(TInt data0, TInt data1){
			return [data0, data1, this](TInt ti){
				TInt ndatatouse = data1 - data0;
				TInt x0 = data0 + (ti*ndatatouse)/this->nthreads;
				TInt x1 = data0 + ((ti+1)*ndatatouse)/this->nthreads;	
				arrutilv2::set_S_H(x1-x0, this->dimension,  this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_sums(), this->get_counts(), this->work_mutex);
			};
		}


		std::vector<std::function<void(TInt)> > base_makeset_C_C_l22s_inds0_mati(TInt ind0, TInt ind1){
			
			TInt ndatatouse = ind1 - ind0;
			
			std::vector<std::function<void(TInt)> > tasks;
			auto set_C_Cl22s_from_inds0etc = [this](TInt ti){
				TInt c0 = (ti*this->ncentroids)/this->nthreads;
				TInt c1 = ((ti+1)*this->ncentroids)/this->nthreads;
				arrutilv2::copyatuniqueindices(this->ndata, c1 - c0, this->dimension, this->data, this->get_C() + c0*this->dimension, this->inds0.get() + c0);
				arrutilv2::set_rl22s(c1 - c0, this->dimension, this->get_C() + c0*this->dimension, this->get_C_l22s() + c0);
			};
		
			if (this->C_init != nullptr && this->initialisation_method.compare("from_C") == 0){
				this->inds0.reset(nullptr);
				
				std::memcpy(this->get_C(), this->C_init, this->ncentroids*this->dimension*sizeof(TFloat));
				arrutilv2::set_rl22s(this->ncentroids, this->dimension, this->get_C(), this->get_C_l22s());
			}
			
			else if (this->data_indices_init_from != nullptr && this->initialisation_method.compare("from_indices") == 0){
				this->inds0 = arrutilv2::copy_ptrarr_to_uptrarr(this->ncentroids, this->data_indices_init_from);
				tasks.push_back(set_C_Cl22s_from_inds0etc);
			}
			
			else if (this->initialisation_method.compare("uniform") == 0){	
				auto C_inds0 = kmeans::initialise2::get_initialisation_indices(this->ncentroids, ndatatouse, this->dimension, this->data + ind0*this->dimension);
				this->C = std::move(std::get<0>(C_inds0));
				arrutilv2::set_rl22s(this->ncentroids, this->dimension, this->get_C(), this->get_C_l22s());
				this->inds0.reset(new TInt [this->ncentroids]);
				auto indiceschosen = std::get<1>(C_inds0).data();
				for (TInt ci = 0; ci < this->ncentroids; ++ ci){
					this->inds0[ci] = ind0 + indiceschosen[ci];
				}
			}
			
			else if (this->initialisation_method.compare("kmeans++") == 0){
				
				auto C_Cl22s_ind0s_mse0 = kmeans::initialise2::get_kmeanspp_initialisation(ndatatouse, this->dimension, this->data + ind0*this->dimension, this->get_data_l22s() + ind0, this->ncentroids);
				this->C = std::move(std::get<0>(C_Cl22s_ind0s_mse0));
				this->C_l22s = std::move(std::get<1>(C_Cl22s_ind0s_mse0));
				this->inds0 = std::move(std::get<2>(C_Cl22s_ind0s_mse0));
				for (TInt ci = 0; ci < this->ncentroids; ++ci){
					this->inds0[ci] += ind0;
				}
				this->mse = std::get<3>(C_Cl22s_ind0s_mse0);
			}
				
			
			else {
				throw std::runtime_error("expected initialisation method in {from_C, from_indices, uniform, kmeans++}. unrecognised initialisation scheme, bailing");
			}
			
			return tasks;
		}
		
	
		virtual void set_initialisation_tasks() = 0;

	
		/*can be used if the principal label, bound etc updater is parameterised as
		TInt this->ncentroids, TInt this->dimension, TFloat * const S, TInt * const H, TInt & this->nchanges, TInt & ndcalcs, TInt this->ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const  C_l22s, otherargs... */
		template <typename Function, typename... Args>
		void base_pll_principal_X(TInt data_start, TInt data_end, const Function & X_updater, TInt ti, Args&&... args){
			TInt ndata = data_end - data_start;
			TInt x0 = data_start + (ti*ndata)/this->nthreads;
			TInt x1 = data_start + ((ti+1)*ndata)/this->nthreads;		
			arrutilv2::pll_update_L_etc(X_updater, 
			this->ncentroids, this->dimension, this->get_sums(), this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_counts(), this->get_dcounts() + ti*this->ncentroids, this->nchanges, this->ndcalcs_X, this->work_mutex, x1-x0, this->data +x0*this->dimension, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), 
			std::forward<Args>(args)...);
		}


		//naive implementation..
		virtual TFloat get_validation_mse(){
			if (this->nvaldata == 0){
				throw std::logic_error("request to compute validataion mse, but nvaldata is 0, bailing");
			}
			
			//sophisticated : 
			
			//get_sse_batchwise(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt & ndcalcs)
		
			TInt local_ndcalcs = 0;	
			TFloat sse = arrutilv2::get_sse_batchwise<TInt, TFloat>(this->nvaldata, this->maxpermultiplyblock_for_validation, this->dimension, this->valdata, this->ncentroids, this->get_C(), this->get_valdata_l22s(), this->get_C_l22s(), local_ndcalcs);
			this->ndcalcs_notX += local_ndcalcs;

			
			return sse/static_cast<TFloat> (this -> nvaldata);
			
								
		}
	
		//this is not exactly what we want for minibatch kmeans, but I will leave it here for now
		virtual void set_n_empty_clusters(){
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
		
		TFloat * const get_dsums(){
			return dsums.get();
		}
		
		TInt * const get_dcounts(){
			return dcounts.get();
		}
	
		virtual void verbose_write_additional(){
			cluster::BaseCluster<TInt, TFloat>::verbose_write_additional();
			
			if (this->round == 0){
				this->verbose_file << "datapoint:\n";
				for (TInt d = 0; d < this->dimension; ++d){
					this->verbose_file << data[d] << "\t";
				}
				this->verbose_file << "\n\n";
			}
			this->verbose_file << "\n********************************************";
			this->verbose_file << "\nthis->round:\n" << this->round << "\n"; 
			this->verbose_file << "\this->ncentroids:\n";
			for (TInt ci = 0; ci < this->ncentroids; ++ci){
				for (TInt d = 0; d < this->dimension; ++d){
					this->verbose_file << C[ci*this->dimension + d] << "\t";
				}
				this->verbose_file << "\n";
			}
			this->verbose_file << "\nlabel:\n" << this->L[0] << "\n";
		}
			

	public:

		virtual ~BaseKmeans(){}

		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const std::string & initialisation_method, const TFloat * const C_init, const TInt * const data_indices_init_from, bool setseed, TInt seed, TFloat maxtime, TInt maxrounds, const std::string & verbose_filename, TInt nvaldata, const TFloat * const valdata, TInt valperiod): 
		
		cluster::BaseCluster<TInt, TFloat> (nthreads, ndata, dimension, ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, data_indices_init_from, setseed, seed, maxtime, maxrounds, verbose_filename, valperiod, nvaldata), data(data), C_init(C_init),  data_l22s(new TFloat [ndata]), C(new TFloat [ncentroids*dimension]), C_l22s(new TFloat [ncentroids]), sums(new TFloat [ncentroids*dimension]), counts(new TInt [ncentroids]), dsums (new TFloat [nthreads*ncentroids*dimension]), dcounts (new TInt [nthreads*ncentroids]), valdata(valdata), valdata_l22s (new TFloat [nvaldata] ){

			this->checkinitialisationvalidity(this->C_init != nullptr);
			std::fill_n(this->get_sums(), this->ncentroids*this->dimension, 0);
			std::fill_n(this->get_counts(), this->ncentroids, 0);
			data_l22s = arrutilv2::get_rl22s(this->ndata, this->dimension, this->data);
			
			this->maxpermultiplyblock_for_validation = std::max(static_cast<TInt> (1),
			static_cast<TInt> ((this->getndata() * this->getdimension())/(2 * this->getncentroids())));


			if (nvaldata > 0){
				valdata_l22s = arrutilv2::get_rl22s(this->nvaldata, this->dimension, this->valdata);
			}
		}
		
		
		
		/* C++11 constructor delegation */
		
		/* Base constructor for no validation */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const std::string & initialisation_method, const TFloat * const C_init, const TInt * const data_indices_init_from, bool setseed, TInt seed, TFloat maxtime, TInt maxrounds, const std::string & verbose_filename):
		BaseKmeans<TInt, TFloat>(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, C_init, data_indices_init_from, setseed, seed, maxtime, maxrounds, verbose_filename, static_cast<TInt> (0), static_cast<const TFloat *> (nullptr), static_cast<TInt> (0) ) {
			}

		
		/* will assume no validation for all the following constructors */
		
		/* from C, no seeding, file_verbosity != 2 */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const TFloat * const C_init, TFloat maxtime, TInt maxrounds) : //TInt seed, 
		BaseKmeans<TInt, TFloat>(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, file_verbosity, file, "from_C", C_init, static_cast <const TInt*> (nullptr), false, static_cast<TInt> (0), maxtime, maxrounds, "") {}
		
		/* from C, with seeding, file_verbosity != 2 (might be some random process in inherited class, even if initialision from C not random) */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const TFloat * const C_init, TInt seed, TFloat maxtime, TInt maxrounds) :  
		BaseKmeans(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, file_verbosity, file, "from_C", C_init, static_cast <const TInt*> (nullptr), true, static_cast<TInt> (0), maxtime, maxrounds, "") {}
		
		
		/* from indices, with seeding, file_verbosity != 2 */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const TInt * const data_indices_init_from, TInt seed, TFloat maxtime, TInt maxrounds) : 
		BaseKmeans(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, file_verbosity, file, "from_indices", static_cast <const TFloat*> (nullptr) , data_indices_init_from, true, seed, maxtime, maxrounds, "") {}

		
		/* from indices, no seeding, file_verbosity != 2 */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const TInt * const data_indices_init_from, TFloat maxtime, TInt maxrounds) : 
		BaseKmeans(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, file_verbosity, file, "from_indices", static_cast <const TFloat*> (nullptr) , data_indices_init_from, false, static_cast<TInt>(0), maxtime, maxrounds, "") {}

		/* initialisation method {uniform, kmeans++ in the future}, with seeding, file_verbosity != 2 */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const std::string & initialisation_method, TInt seed, TFloat maxtime, TInt maxrounds) : 
		BaseKmeans(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, static_cast <const TFloat*> (nullptr), static_cast <const TInt*> (nullptr), true, seed, maxtime, maxrounds, "") {}
		
		
		/* version usable by the hulk with option for nthreads and extreme verbosity, no seeding */
		BaseKmeans(int cout_verbosity, TInt ndata, TInt dimension, const TFloat *  const data, TInt ncentroids, TFloat * const C, std::ofstream & file, TFloat maxtime, TInt nthreads, const std::string & verbose_filename) : 
		BaseKmeans<TInt, TFloat>(nthreads, ndata, dimension, data, ncentroids, cout_verbosity, 2, file, "from_C", C, static_cast <const TInt*> (nullptr), false, static_cast<TInt> (0), maxtime, std::numeric_limits<TInt>::max(), verbose_filename) {}
		
		
		/*maxtime and maxround as large as possible, (version used in clustering centroids, see pllinityinyang) */
		BaseKmeans(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, 
		int verbosity, std::ofstream & file, const TInt * const data_indices_init_from): BaseKmeans<TInt, TFloat>(nthreads, ndata, dimension, data, ncentroids, verbosity, 0, file, "from_indices", static_cast <const TFloat*> (nullptr), data_indices_init_from, false, static_cast<TInt> (0), std::numeric_limits<TFloat>::max(), std::numeric_limits<TInt>::max(), "") {}






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
		
		
		
		
		/* return : data_rearranged, grouparts, groupsizes (inspired by .../clustering/kmeansfunctionwrapper.h)  */
		std::tuple< std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> > 
		/* where minclustersize is the minimum number of data per cluster, enforced by random reassignemnts after clustering. Useful for yinyang, TODO: incorporate equalkmeans code here */
		get_contig_by_cluster_3(TInt minclustersize = 0){
			
			/* do the clustering */
			auto C_L_ind0 = this->get();
			auto L = std::move(std::get<1>(C_L_ind0));
			
			/* get counts and make changes if nec to guarantee that minclustersize recognised */
			std::unique_ptr< TInt []> groupsizes = arrutilv2::gethistogram(this->ndata, this->ncentroids, L.get());
			arrutilv2::make_balanced(minclustersize, this->ndata, L.get(), this->ncentroids, groupsizes.get());
			/* get parts is simply obtained from groupsizes */
			std::unique_ptr< TInt []> groupparts (new TInt [this->ncentroids + 1] );
			groupparts[0] = 0;
			for (TInt pi = 0; pi < this->ncentroids; ++pi){
				groupparts[pi + 1] = groupparts[pi] + groupsizes[pi];
			}
			
			/* do the rearrangement */			
			std::vector<TInt> currentcounts (this->ncentroids, 0);
			std::unique_ptr< TFloat []> data_rearranged (new TFloat [this->ndata*this->dimension] );
			for (TInt i = 0; i < this->ndata; ++i){
				std::memcpy(data_rearranged.get() + (groupparts[L[i]] + currentcounts[L[i]])*this->dimension, data + i*this->dimension, sizeof(TFloat)*this->dimension);
				++ currentcounts[L[i]];
			}
			return std::make_tuple(std::move(data_rearranged), std::move(groupparts), std::move(groupsizes));
		}
		
		
		
		/* return : version of above, but additionally returns data_rearranged_l22s (paramater 1) */
		std::tuple< std::unique_ptr<TFloat []>, std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> > 
		get_contig_by_cluster_4(TInt minclustersize = 0){
			auto datr_gp_gs = this->get_contig_by_cluster_3(minclustersize);
			std::unique_ptr<TFloat []> datr_l22s (new TFloat [this->ndata]);
			arrutilv2::set_rl22s(this->ndata, this->dimension, std::get<0>(datr_gp_gs).get(), datr_l22s.get()); 
			return std::make_tuple(std::move(std::get<0>(datr_gp_gs)), std::move(datr_l22s), std::move(std::get<1>(datr_gp_gs)), std::move(std::get<2>(datr_gp_gs)));
		}

		const TFloat * const get_data_l22s(){
			return data_l22s.get();
		}

		const TFloat * const getdata(){
			return data;
		}
		
		const TFloat * const get_C_init(){
			return C_init;
		}

		//Return an estimate of memory requirement of this base class
		virtual TInt get_approximate_memory_requirement(){
			return sizeof(TFloat)*(
			this->ncentroids*this->dimension + // C
			this->ncentroids +// C_l22s
			this->ncentroids*this->dimension + //sums
			this->ncentroids*this->dimension*this->nthreads + //dsums
			this->ndata*this->dimension + //data (not really part of object but hey)
			this->nvaldata + //valdata_l22s
			this->ndata) + //data_l22s

			
			sizeof(TInt)*(
			this->ndata + // L
			this->ncentroids + //counts
			this->ncentroids*this->nthreads);	//dcounts	
			
		}
		
};

}


#endif


