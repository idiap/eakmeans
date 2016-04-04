#ifndef PLL_BASEKMEANSSS_H
#define PLL_BASEKMEANSSS_H

#include "basedensecentroidkmeans.h"

#include "arrutilv2l3.h"
#include "arrutilv2mse.h"
#include "arrutilv2discrete.h"

#include "sample.h"
#include "initialise2.h"





namespace kmeans{

template <typename TInt, typename TFloat>
class BaseKmeans : public kmeans::BaseDenseCentroidKmeans<TInt, TFloat> { 

	private:
		virtual std::function<void(TInt)> update_L_S_H_ati(){
			throw std::logic_error("Call to unimplemented version update_L_S_H_ati in BaseKmeans, unexpected. Perhaps the call is from a derived class which should not update these variables simultaneously? Bailing.");
			return [](TInt){};
		}
		
	
	
			/* what to do with data samples which were used in the previous round (don't need any initialisation) used by grow batch algorithms*/
		virtual void update_already_used(TInt x0, TInt x1, TInt ti){
			throw std::logic_error("Call to update_already_used (in basedensecentroidkmeans) is unexpected. This function is only for brow batch algorithms");
		};
		
		/* what to do with data samples which were not used in the previous round (need some type of initialisation) used by grow batch algorithms */
		virtual void update_unused(TInt x0, TInt x1, TInt ti){
			throw std::logic_error("Call to update_ununsed (in basedensecentroidkmeans) is unexpected. This function is only for brow batch algorithms");
		};
	
		
	protected:
	
		const TFloat * data;
		std::unique_ptr<TFloat []> dsums;
		std::unique_ptr<TInt []> dcounts;
		const TFloat * valdata;
		TInt maxpermultiplyblock_for_validation;

		virtual void set_L(TInt x0, TInt x1) override final {
			TInt local_ndcalcs = 0;			
			arrutilv2::set_rargmins(x1 - x0, this->dimension, this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, local_ndcalcs);
			this->ndcalcs_X += local_ndcalcs;

		} 
		
		virtual std::function<void(TInt)> base_set_S_H_ati(TInt data0, TInt data1) override final{
			return [data0, data1, this](TInt ti){
				TInt ndatatouse = data1 - data0;
				TInt x0 = data0 + (ti*ndatatouse)/this->nthreads;
				TInt x1 = data0 + ((ti+1)*ndatatouse)/this->nthreads;	
				this->set_S_H(x0, x1, ti);
			};
		}

		virtual void set_S_H(TInt x0, TInt x1, TInt ti){
			arrutilv2::set_S_H(x1-x0, this->dimension,  this->data + x0*this->dimension, this->ncentroids, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_dcounts() + ti*this->ncentroids,  this->get_sums(), this->get_counts(), this->work_mutex);
		}		

		virtual std::function<void(TInt)> set_C_Cl22s_from_inds0etc_ati(){
		
			return [this](TInt ti){
				
				TInt c0 = (ti*this->ncentroids)/this->nthreads;
				TInt c1 = ((ti+1)*this->ncentroids)/this->nthreads;
				arrutilv2::copyatuniqueindices(this->ndata, c1 - c0, this->dimension, this->data, this->get_C() + c0*this->dimension, this->inds0.get() + c0);
				arrutilv2::set_rl22s(c1 - c0, this->dimension, this->get_C() + c0*this->dimension, this->get_C_l22s() + c0);
			};
			
		}
		
		virtual std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt [] > > get_C_inds0_uniform(TInt ind0, TInt ind1){


			TInt ndatatouse = ind1 - ind0;
			auto C_inds0 = kmeans::initialise2::get_initialisation_indices(this->ncentroids, ndatatouse, this->dimension, this->data + ind0*this->dimension);
			
			std::unique_ptr < TInt []> inds0 (new TInt [this->ncentroids]);
			auto indiceschosen = std::get<1>(C_inds0).data();
			for (TInt ci = 0; ci < this->ncentroids; ++ ci){
				inds0[ci] = ind0 + indiceschosen[ci];
			}
			
			return std::make_tuple (std::move(std::get<0>(C_inds0)), std::move(inds0));
		}
		
		virtual void do_kmeanspp_initialisation(TInt ind0, TInt ind1){
			

			TInt ndatatouse = ind1 - ind0;
				
			auto C_Cl22s_ind0s_mse0 = kmeans::initialise2::get_kmeanspp_initialisation(ndatatouse, this->dimension, this->data + ind0*this->dimension, this->get_data_l22s() + ind0, this->ncentroids);
			
			
			std::memcpy(this->C.get(), std::get<0>(C_Cl22s_ind0s_mse0).get(), sizeof(TFloat)*this->ncentroids*this->dimension);
			//this->C = std::move(std::get<0>(C_Cl22s_ind0s_mse0));
			

			std::memcpy(this->C_l22s.get(), std::get<1>(C_Cl22s_ind0s_mse0).get(), sizeof(TFloat)*this->ncentroids);
			//this->C_l22s = std::move(std::get<1>(C_Cl22s_ind0s_mse0));
			
			//this->inds0 = std::move(std::get<2>(C_Cl22s_ind0s_mse0));			
			std::memcpy(this->inds0.get(), std::get<2>(C_Cl22s_ind0s_mse0).get(), sizeof(TInt)*this->ncentroids);

			for (TInt ci = 0; ci < this->ncentroids; ++ci){
				this->inds0[ci] += ind0;
			}
			this->mse = std::get<3>(C_Cl22s_ind0s_mse0);
		}
		
		
		virtual TFloat getmeanl22at() override final{
			return arrutilv2::getmeanl22at(this->ncentroids, this->dimension, this->get_C(), this->ndata, this->data, this->get_L(), this->get_C_l22s(), this->get_data_l22s());
		}	

		
		//used by certain descendants. I believe that this eventually uses blas matrix multiplication (if withblas). 
		virtual void set_upper_lowers_L(TInt x0, TInt x1) override final{ /* from basedensecentroidkmeans */
			TInt local_ndcalcs = 0;
			arrutilv2::set_rrl2ss_argminmins(x1 - x0, this->getdimension(), this->getdata() + x0*this->getdimension(), this->getncentroids(), this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->elkan_lowers_base.get() + x0*this->getncentroids(), this->get_L() + x0, this->elkan_upper_base.get() + x0, local_ndcalcs);
			this->ndcalcs_X += local_ndcalcs;
		}
		


		/* used by dense grow batch algorithms. Works fine if only one task per batch.  */
		std::vector<std::function<void(TInt)> > bgbmse_update_L_dn_etc_S_H_batch_switch_mati(const growbatchapp::GBApp<TInt, TFloat> & gba){
		std::vector<std::function<void(TInt)> > tasks = {};
		tasks.emplace_back (
			//update L and dn of data used in previous round and update S and H accordingly. 				
			[this, &gba](TInt ti){
				TInt x0 = (ti*gba.ndata_active_previous)/this->nthreads;
				TInt x1 = ((ti+1)*gba.ndata_active_previous)/this->nthreads;		
				this->update_already_used(x0, x1, ti);
			}
		);
			
		tasks.emplace_back(
			[this, &gba](TInt ti){
			//if there is additional data in this round, it should be labelled and S and H should be updated accordingly. 
				if (gba.ndata_active != gba.ndata_active_previous){
					TInt ndata_tail = gba.ndata_active - gba.ndata_active_previous;
					TInt x0 =  gba.ndata_active_previous + (ti*ndata_tail)/this->nthreads;
					TInt x1 =  gba.ndata_active_previous + ((ti + 1)*ndata_tail)/this->nthreads;
					this->update_unused(x0, x1, ti);
				} 
			}
		);
		return tasks;
	}


		
	
		virtual void set_initialisation_tasks() = 0;

	
		/* Can be used if the principal label, bound etc updater is parameterised as
		TInt this->ncentroids, TInt this->dimension, TFloat * const S, TInt * const H, TInt & this->nchanges, TInt & ndcalcs, TInt this->ndata, const TFloat * const data, const TFloat * const C, const TFloat * const data_l22s,  const TFloat * const  C_l22s, otherargs... */
		//The use of data_start and data_end is confusing.
		template <typename Function, typename... Args>
		void base_pll_principal_X(TInt data_start, TInt data_end, const Function & X_updater, TInt ti, Args&&... args){ // x
			TInt ndata = data_end - data_start;
			TInt x0 = data_start + (ti*ndata)/this->nthreads;
			TInt x1 = data_start + ((ti+1)*ndata)/this->nthreads;		
			
			arrutilv2::pll_update_L_etc(
			//The compulsory parameters to pll_update_L_etc,
			X_updater, 
			this->ncentroids, this->dimension, this->get_sums(), this->get_dsums() + ti*this->dimension*this->ncentroids, this->get_counts(), this->get_dcounts() + ti*this->ncentroids, this->nchanges, this->ndcalcs_X, this->work_mutex, 
			//The additional parameters to pll_update_L_etc:
			x1-x0, this->data +x0*this->dimension, this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), 
			std::forward<Args>(args)...);
		}
		
		


		//naive implementation..
		virtual TFloat get_validation_mse(){
			if (this->nvaldata == 0){
				throw std::logic_error("request to compute validataion mse, but nvaldata is 0, bailing");
			}
		
			TInt local_ndcalcs = 0;	
			TFloat sse = arrutilv2::get_sse_batchwise<TInt, TFloat>(this->nvaldata, this->maxpermultiplyblock_for_validation, this->dimension, this->valdata, this->ncentroids, this->get_C(), this->get_valdata_l22s(), this->get_C_l22s(), local_ndcalcs);
			this->ndcalcs_notX += local_ndcalcs;

			
			return sse/static_cast<TFloat> (this -> nvaldata);
			
								
		}
		
		
		virtual std::unique_ptr<TFloat []> getnew_data_l22s() override final{
			
			return arrutilv2::get_rl22s(this->ndata, this->dimension, this->data);
		}
		
		virtual std::unique_ptr<TFloat []> getnew_valdata_l22s() override{
			return arrutilv2::get_rl22s(this->nvaldata, this->dimension, this->valdata);
		}


	

		virtual void set_X_tasks() = 0;
		
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
					this->verbose_file << this->C[ci*this->dimension + d] << "\t";
				}
				this->verbose_file << "\n";
			}
			this->verbose_file << "\nlabel:\n" << this->L[0] << "\n";
		}





		virtual TFloat get_clean_mse(bool on_validation) override final{
			//set data to use depending on on_validation
			auto cmse_data = this->data;
			auto cmse_data_l22s = this->get_data_l22s();
			auto cmse_ndata = this->ndata;
		
			if (on_validation == true){
				cmse_data = this->valdata;
				cmse_data_l22s = this->valdata_l22s.get();
				cmse_ndata = this->nvaldata;
			}
			
			TFloat sse = 0;
			TInt argmin; 
			TFloat min;
			std::unique_ptr<TFloat []> l2s ( new TFloat [this->ncentroids]);
			for (TInt i = 0; i < cmse_ndata; ++i){
				arrutilv2::set_rl2s_argminmin(this->dimension, cmse_data + i*this->dimension, this->ncentroids, this->get_C(), cmse_data_l22s[i], this->get_C_l22s(), l2s.get(), argmin, min);
				sse += min*min;
			}
			return sse / static_cast<TInt> (this->ndata);
		}

			

	public:

		virtual ~BaseKmeans(){}


	
		
		BaseKmeans(
		TInt nthreads, //1
		TInt ndata, //2
		TInt dimension, //3 
		const TFloat * const data, //4 *** not forwarded to parent
		TInt ncentroids, //5
		int cout_verbosity, //6
		int file_verbosity, //7
		std::ofstream & file, //8
		const std::string & initialisation_method, //9
		const TFloat * const C_init, //10
		const TInt * const data_indices_init_from, //11
		bool setseed, //12
		TInt seed, //13
		TFloat maxtime, //14
		TInt maxrounds, //15
		const std::string & verbose_filename, //16
		TInt nvaldata, //17
		const TFloat * const valdata, //18 *** not forwarded to parent
		TInt valperiod, //19
		const std::string & cmsewritefn, //20
		TInt cmserate,	//21
		TFloat gbphi //22
		):
		
		kmeans::BaseDenseCentroidKmeans<TInt, TFloat> (
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
		gbphi,
		C_init //20
		), 
			
		data(data), dsums (new TFloat [nthreads*ncentroids*dimension]), dcounts (new TInt [nthreads*ncentroids]), valdata(valdata){

			this->maxpermultiplyblock_for_validation = std::max(static_cast<TInt> (1),
			static_cast<TInt> ((this->getndata() * this->getdimension())/(2 * this->getncentroids())));

	}


		/* Using C++11 constructor delegation, the following is the version used in clustering centroids, see pllinityinyang)
		 * maxtime and maxround as large as possible, etc etc */
		
		BaseKmeans(
		TInt nthreads, //1
		TInt ndata, //2
		TInt dimension, //3 
		const TFloat * const data, //4
		TInt ncentroids, //5
		int verbosity, //6
		std::ofstream & file, //7 
		const TInt * const data_indices_init_from  //8
		):
		 
		BaseKmeans<TInt, TFloat>(
		nthreads, //1
		ndata, //2
		dimension, //3 
		data, //4
		ncentroids, //5 
		verbosity, //6
		0, //7
		file, //8
		"from_indices", //9
		static_cast <const TFloat*> (nullptr), //10
		data_indices_init_from, //11
		false, //12
		static_cast<TInt> (0), //13
		std::numeric_limits<TFloat>::max(), //14
		std::numeric_limits<TInt>::max(), //15
		"", //16
		static_cast<TInt> (0), //17
		static_cast<const TFloat *> (nullptr), //18 
		static_cast<TInt> (0), //19
		"",
		static_cast<TInt> (0),
		static_cast<TFloat> (0.0)
		) {}





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


		const TFloat * const getdata(){
			return data;
		}
		
		const TFloat * const getvaldata(){
			return valdata;
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

