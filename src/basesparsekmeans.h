#ifndef PLL_BASEKDENSECENTROIDKMEANS_H
#define PLL_BASEKDENSECENTROIDKMEANS_H

#include "basedensecentroidkmeans.h"
#include "arrutilv2l3.h"
#include "arrutilv2copy.h"
#include "arrutilv2mse.h"
#include "arrutilv2discrete.h"
#include "sample.h"
#include "stringutilfile.h"
#include "sparsedatasets.h"
#include "sparseutil.h"
#include "sparseinitialise.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseSparseKmeans : public kmeans::BaseDenseCentroidKmeans<TInt, TFloat> {

	private:
		std::unique_ptr<sparse::SparseData<TInt, TFloat> > private_data; //used when data read from file, otherwise empty
		std::unique_ptr<sparse::SparseData<TInt, TFloat> > private_valdata; //used when data read from file, otherwise empty

	protected:
	
	
		virtual TFloat get_clean_mse(bool on_validation) override final{
			
			
						
			//set data to use depending on on_validation
			auto cmse_ptrdata = this->ptrdata;
			auto cmse_data_l22s = this->get_data_l22s();
			auto cmse_ndata = this->ndata;
		
			if (on_validation == true){
				cmse_ptrdata = this->ptrvaldata;
				cmse_data_l22s = this->valdata_l22s.get();
				cmse_ndata = this->nvaldata;
			}
			
			
			TFloat sse = 0;
			TInt argmin; 
			TFloat min;
			std::unique_ptr<TFloat []> l2s ( new TFloat [this->ncentroids]);
			for (TInt i = 0; i < cmse_ndata; ++i){
				sparse::set_argminmin_rl2s(
				cmse_ptrdata->starts[i+1] - cmse_ptrdata->starts[i], 
				cmse_ptrdata->indices.data() + cmse_ptrdata->starts[i], 
				cmse_ptrdata->values.data() + cmse_ptrdata->starts[i], 
				this->dimension, 
				this->ncentroids, 
				this->get_C(), 
				cmse_data_l22s[i], 
				this->get_C_l22s(), 
				argmin,
				min, 
				l2s.get());
				
				sse += min*min;
			}
			
						 
			return sse / static_cast<TInt> (this->ndata);
				
		
		}


		virtual TFloat getmeanl22at() override final{
			return sparse::getmeanl22at(*this->ptrdata, this->get_C(), this->get_L(), this->get_data_l22s(), this->get_C_l22s());			
		}	
		
		virtual void verbose_write_additional(){
			throw std::runtime_error("verbose_write_additional needs implementing in basesparsekmeans");
		}
		
		virtual void set_mse() = 0;
		
		virtual void set_L(TInt x0, TInt x1) override final{
			
			for (TInt i = x0; i < x1; ++i){
				sparse::set_label(this->ptrdata->starts[i+1] - this->ptrdata->starts[i], 
				this->ptrdata->indices.data() + this->ptrdata->starts[i], 
				this->ptrdata->values.data() + this->ptrdata->starts[i], this->ptrdata->dimension, this->ncentroids, 
				this->get_C(), //this->data_l22s[i], 
				this->get_C_l22s(), this->L[i]);
			}
		}
			
	
		const sparse::SparseData<TInt, TFloat> * ptrdata;
		const sparse::SparseData<TInt, TFloat> * ptrvaldata;
		
		//will be of length nthreads, each vector containing (index, old, new) where index is index of data whose label changed
		std::vector<std::vector<std::tuple<TInt, TInt, TInt > > > where_label_changes; 
		


		

		void set_S_H(TInt data0, TInt data1){
			sparse::todense::set_S_H(
			*this->ptrdata, data0, data1, this->ncentroids, this->get_L(),  this->get_sums(),  this->get_counts()			
			);
		}
		
		
		
		//A hack as no pllsation as suggested by ati suffix
		virtual std::function<void(TInt)> base_set_S_H_ati(TInt data0, TInt data1) override final{
			return [data0, data1, this](TInt ti){
				if (ti == 0){
					this->set_S_H(data0, data1);
				}
			};
		}



		virtual std::function<void(TInt)> set_C_Cl22s_from_inds0etc_ati() override {
			return [this](TInt ti){				
				TInt c0 = (ti*this->ncentroids)/this->nthreads;
				TInt c1 = ((ti+1)*this->ncentroids)/this->nthreads;
				
				//for (TInt ci = 0; ci < this->ncentroids; ++ ci){
					//std::cout << this->inds0[0] << " " << std::flush;
				//}
				std::cout << std::endl;
				sparse::todense::copyatindices(c1 - c0, *this->ptrdata, this->get_C() + c0*this->dimension, this->inds0.get() + c0);
				arrutilv2::set_rl22s(c1 - c0, this->dimension, this->get_C() + c0*this->dimension, this->get_C_l22s() + c0);

			};
		}
		
		virtual std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt [] > > get_C_inds0_uniform(TInt ind0, TInt ind1){
			auto C_inds0 = kmeans::sparseinit::get_initialisation_indices(this->ncentroids, *this->ptrdata, ind0, ind1);
			auto inds0 = arrutilv2::get_with_offset(this->ncentroids, std::get<1>(C_inds0).data(), ind0);			
			return std::make_tuple ( std::move(std::get<0>(C_inds0)), std::move(inds0) );
		}
		
		virtual void do_kmeanspp_initialisation(TInt ind0, TInt ind1){
			
			auto C_Cl22s_ind0s_mse0 = kmeans::sparseinit::get_kmeanspp_initialisation(this->ncentroids, *this->ptrdata, ind0, ind1);
			this->C = std::move(std::get<0>(C_Cl22s_ind0s_mse0));
			this->C_l22s = std::move(std::get<1>(C_Cl22s_ind0s_mse0));
			this->inds0 = std::move(std::get<2>(C_Cl22s_ind0s_mse0));
			
			for (TInt ci = 0; ci < this->ncentroids; ++ci){
				this->inds0[ci] += ind0;
			}
			this->mse = std::get<3>(C_Cl22s_ind0s_mse0);
			
		}
		
		
	
		virtual void set_initialisation_tasks() = 0;



		//untested
		virtual TFloat get_validation_mse(){
			

			if (this->nvaldata == 0){
				throw std::logic_error("request to compute validataion mse, but nvaldata is 0, bailing");
			}
			
			
			TFloat sse = 0;
			TInt label = 0;
			
			for (TInt i = 0; i < this->ndata; ++i){
				
				//TODO: the following duplicates computation (find label, then compute distance to label)
				
				sparse::set_label(this->ptrvaldata->starts[i+1] - this->ptrvaldata->starts[i], 
				this->ptrvaldata->indices.data() + this->ptrvaldata->starts[i], 
				this->ptrvaldata->values.data() + this->ptrvaldata->starts[i],
				this->dimension, this->ncentroids, this->get_C(), //this->valdata_l22s[i], 
				this->get_C_l22s(), label);
					
				
				sse += 
				this->valdata_l22s[i] + this->C_l22s[label]				
				-2.*sparse::get_inner(this->ptrvaldata->starts[i+1] - this->ptrvaldata->starts[i], 
				this->ptrvaldata->indices.data() + this->ptrvaldata->starts[i], 
				this->ptrvaldata->values.data() + this->ptrvaldata->starts[i],
				this->get_C() + this->dimension*label);
			
			}
			return sse / static_cast<TFloat> (this->ndata);
			

		}
		
		
		virtual std::unique_ptr<TFloat []> getnew_data_l22s() override{
			return sparse::get_rl22s(*this->ptrdata);
		}
		
		virtual std::unique_ptr<TFloat []> getnew_valdata_l22s() override{
			return sparse::get_rl22s(*this->ptrvaldata);
		}


		virtual void set_X_tasks() = 0;
		
	
		
		
	public:

		virtual ~BaseSparseKmeans(){}

		//initialise from file.
		BaseSparseKmeans(
		TInt nthreads, //1
		TInt ndata, //2
		TInt dimension, //3 
		const std::string & datafin, //4 ***
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
		const std::string & valdatafin, //18  ***
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
		gbphi, //19
		C_init //20
		), 
		
		where_label_changes (this->nthreads, std::vector<std::tuple<TInt, TInt, TInt>> {} )
		
		
		{
			if (datafin.compare("") != 0){
				bool withheader = stringutilfile::file_has_2int_header(datafin);
				this->private_data.reset(new sparse::SparseData<TInt, TFloat>(datafin, withheader));
				this->ptrdata = private_data.get();
			}
			
			if (valdatafin.compare("") != 0){
				bool withheader = stringutilfile::file_has_2int_header(valdatafin);
				this->private_valdata.reset(new sparse::SparseData<TInt, TFloat>(valdatafin, withheader));
				this->ptrvaldata = private_valdata.get();
			}
		}
		
		//initialise from SparseDatas. Same as initializer from files, but strings of filenames are replaced with sparsedatas.
		BaseSparseKmeans(TInt nthreads, TInt ndata, TInt dimension, const sparse::SparseData<TInt, TFloat> * ptrdata, TInt ncentroids, int cout_verbosity, int file_verbosity, std::ofstream & file, const std::string & initialisation_method, const TFloat * const C_init, const TInt * const data_indices_init_from, bool setseed, TInt seed, TFloat maxtime, TInt maxrounds, const std::string & verbose_filename, TInt nvaldata, const sparse::SparseData<TInt, TFloat> * ptrvaldata, TInt valperiod, const std::string & cmsewritefn, TInt cmserate, TFloat gbphi): 
		
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
		gbphi, //19
		C_init //20
		), 
		
		where_label_changes (this->nthreads, std::vector<std::tuple<TInt, TInt, TInt>> {} )
		{
			this->ptrdata = ptrdata;
			this->ptrvaldata = ptrvaldata;
		}
		
		
		
		

		const sparse::SparseData<TInt, TFloat> * const getptrdata(){
			return ptrdata;
		}
		
		const sparse::SparseData<TInt, TFloat> * const getptrvaldata(){
			return ptrvaldata;
		}
		

		//Return an estimate of memory requirement of this base class
		virtual TInt get_approximate_memory_requirement() {
			//TODO
			return 1;
		}
		
};

}


#endif
