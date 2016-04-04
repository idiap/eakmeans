#include "pllcluster.h"
#include "pllkmeansfuncs_nonvoid.h"
#include "pllkmeansfuncs_void.h"
#include "pllkmeansfuncs.hpp"
#include <iostream>
#include <exception>
#include <iterator>


#include "sparsedatasets.h"

#include "stringutilfile.h"
#include "txtdatasets.h"
	
namespace cluster{
	
	template <typename TFloat>
	void set_data_from_file(const std::string &fn, std::unique_ptr<TFloat []> & data, size_t & ndata, size_t & dimension, bool shuffle){
		bool withheader = stringutilfile::file_has_2int_header(fn);
		std::fstream valdfile(fn, std::ios_base::in);
		if (withheader == false){
			std::cout << "Determining #data and dimension from header less file ... " << std::endl;
			std::ifstream myfile(fn);
			myfile.unsetf(std::ios_base::skipws);
			ndata = std::count(std::istream_iterator<char>(myfile), std::istream_iterator<char>(), '\n');
			std::cout << "#data : " << ndata << std::endl;
			std::ifstream dfile(fn, std::ios_base::in);		
			std::string frag;
			std::string aline;
			std::getline(dfile, aline);
			std::stringstream ss(aline);
			dimension = 0;
			while (ss){
				ss >> frag;
				if (ss){ 
					dimension += 1;
				}
			}
			std::cout << "dimension : " << dimension << std::endl;
		}
		
		else{
			valdfile >> ndata;
			valdfile >> dimension;
		}
		
		data.reset(new TFloat [ndata*dimension]);
		TFloat * raw_dataptr = data.get();
		while (valdfile >> *raw_dataptr){
			++raw_dataptr;
		}
		valdfile.close();
		
		if (shuffle == true){
			std::cout << "Shuffling data..." << std::endl;
			randomutil::sample::inplace_shuffle_by_row(ndata, dimension, data.get()); 
			std::cout << "done" << std::endl;
		}
		
		std::cout << "Leaving set_data_from_file " << std::endl;
	}



			//if (dimdata != dimension){
				//throw std::runtime_error("dimension determined by file and that passed to function do not agree. OLD ERROR MESSAGE: The dimension of the training data and validation data do not seem to have the same dimension. Recall the format of the data files : first line is (ndata dimension) followed by ndata lines with dimension floating point values"); 
			//}
			
	
	std::tuple<std::unique_ptr<float []>, std::unique_ptr<size_t []>, std::unique_ptr<size_t []>, size_t, size_t, float, std::string>
	solveiolessf(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const float * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const float * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, float maxtime, size_t maxrounds, size_t minibatchsize, bool captureverbose){
		//std::string emptystring = "";
		return solveioless<float>(algorithm, nthreads, ndata, dimension, data, ncentroids, cout_verbosity, initialisation_method, C_init, data_indices_init_from, setseed, seed, maxtime, maxrounds, minibatchsize, 0, nullptr, 0, captureverbose);
	}
	
	std::tuple<std::unique_ptr<double []>, std::unique_ptr<size_t []>, std::unique_ptr<size_t []>, size_t, size_t, double, std::string>
	solveiolessd(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const double * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const double * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, double maxtime, size_t maxrounds, size_t minibatchsize, bool captureverbose){
		//std::string emptystring = "";
		return solveioless<double>(algorithm, nthreads, ndata, dimension, data, ncentroids, cout_verbosity, initialisation_method, C_init, data_indices_init_from, setseed, seed, maxtime, maxrounds, minibatchsize, 0, nullptr, 0, captureverbose);
	}
	


	template <typename TFloat>
	void solvewrite(
	const std::string & algorithm, 
	bool issparse, 
	size_t nruns, 
	size_t nthreads, 
	int cout_verbosity, 
	int file_verbosity, 
	const std::string & datainfn,
	const std::string & coutfn, 
	const std::string & loutfn,  
	const std::string & ioutfn, 
	const std::string & soutfn, 
	const std::string & voutfn, 
	const std::string & moutfn, 
	const std::string & moutdir, 
	const std::string & cinf, 
	const std::string & ind0fn, 
	const std::string & init0, 
	bool setseed, 
	size_t seed, 
	size_t ncentroids, 
	size_t maxiter, 
	double maxtime, 
	const std::string & valinfn, 
	size_t valperiod, 
	size_t minibatchsize, 
	std::string & cmsewritefn, 
	size_t cmserate, //27
	TFloat gbphi
	){
		
		std::cout << "Just entered solvewrite, with issparse " << issparse << std::endl;
		
		typedef std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<size_t []>, std::unique_ptr<size_t []>, size_t, size_t, TFloat> tup6;
		
		
		size_t nvaldata = 0;
		
		/* will not be used if issparse */
		std::unique_ptr<TFloat []> valdata;
		std::unique_ptr<TFloat []> data;
		
		/* will not be used if not issparse */
		sparse::SparseData<size_t, TFloat> s_data {};
		sparse::SparseData<size_t, TFloat> s_valdata {};
		
		

		size_t ndata;
		size_t dimension;		
		
		/* load data, not issparse version */
		if (!issparse){
			std::cout << "Non-sparse data being loaded into array... " << std::flush;
			
			bool shuffle = setseed;
			
			set_data_from_file<TFloat>(datainfn, data, ndata, dimension, shuffle);
			
			std::cout << "done." << std::endl;
		}
		
		/* load data, issparse version */
		else{
			
			std::string full_datainfn = datainfn;
			
			/* preliminary test to see if file exists, try to find match if not */
			std::fstream dfile(datainfn, std::ios_base::in);
			if (!dfile.is_open()){
				/* try on fixed path TODO: this old hack may no longer work */
				dfile.open(datasets::sparse_data_dir + datainfn);				
				if (!dfile.is_open()){
					throw std::runtime_error("sparse datafile " + datainfn + " not found, neither as is nor on fixed path  " + datasets::sparse_data_dir + " ( " + datasets::sparse_data_dir + datainfn  + " ) " );
				}
				else{
					full_datainfn = datasets::sparse_data_dir + datainfn;
				}
			}
			dfile.close();
			/* preliminary test complete */
	

			bool withheader = stringutilfile::file_has_2int_header(full_datainfn);
			s_data = sparse::SparseData<size_t, TFloat>(full_datainfn, withheader);	
			ndata = s_data.ndata;
			dimension = s_data.dimension;
		}
		
		
			
		
		std::string initialisation_method = init0;
		std::unique_ptr<TFloat []> C_init;
		std::unique_ptr<size_t []> data_indices_init_from;
		if (cinf.compare("") != 0){
			initialisation_method = "from_C";
			C_init.reset(new TFloat[ncentroids*dimension]);
			std::fstream Cfile(cinf, std::ios_base::in);
			TFloat * raw_C_init_ptr = C_init.get();
			while (Cfile >> *raw_C_init_ptr){
				++raw_C_init_ptr;
			}
			Cfile.close();
		}
		
		else if (ind0fn.compare("") != 0){
			
			initialisation_method = "from_indices";
			data_indices_init_from.reset(new size_t[ncentroids]);
			std::fstream init_file(ind0fn, std::ios::in);
			if (init_file.is_open() == false){
				std::cout << "failed to open " << ind0fn << std::endl; 
				throw std::runtime_error("The file of initialising indices did not open successfully");
			}
			
			std::cout << "init_file " << ind0fn << " is open " << std::endl; 
				
			size_t * rp = data_indices_init_from.get();
			while (init_file >> *rp){
	
				++rp;
			}
			init_file.close();	
		}
		
		else {
			
		}
		

		if (valinfn.compare("") != 0){
			//populate validation data, dense version
			if (!issparse){
				size_t val_dimension;
				set_data_from_file<TFloat>(valinfn, valdata, nvaldata, val_dimension, false);
				
				if (val_dimension != dimension){
					throw std::runtime_error("dimension determined by validation file and data file not in agreememt.OLD ERROR MESSAGE: The dimension of the training data and validation data do not seem to have the same dimension. Recall the format of the data files : first line is (ndata dimension) followed by ndata lines with dimension floating point values"); 
				}
			}
			
			else{
				bool withheader = stringutilfile::file_has_2int_header(valinfn);
				s_valdata = sparse::SparseData<size_t, TFloat>(valinfn, withheader);
				nvaldata = s_valdata.ndata;
				size_t dimension_valdata = s_valdata.dimension;

				if (dimension_valdata != dimension){
					
					std::string rem = "dimension determined from training and test data differ. This is not unexpected: a test datapoint has a non-zero value where no training data in non-zero. At this point in the code (in pllkmeansfuncs.cpp) it suffices to set dimension to max (dimension, dimension_valdata) to fix the problem which may arise of dimension is left as is. I suggest you come into the code at this point and uncomment the line just below! (Alternatively, what about passing in files with dimensions as headers? I think these exist and are ready to use...)";
					
					rem = rem + "\ndim from train : " + std::to_string(dimension) + "\t dim from test : " + std::to_string(dimension_valdata) + "\n";
					
					//throw std::runtime_error(rem);
				}
	
				dimension = std::max(dimension_valdata, dimension);
				s_valdata.dimension = dimension;
				s_data.dimension = dimension;
				 
			}
		}
		
		
		
		//best training mse. 
		TFloat bestmse = std::numeric_limits<TFloat>::max();
		for (size_t ri = 0; ri < nruns; ++ri){
			
			std::ofstream file;
			if (soutfn.compare("") != 0){
				file.open(soutfn);		
			}
			
			else if (moutdir.compare("") != 0){
				file.open(moutdir + "run" + std::to_string(ri) + ".txt");
			}
			
			tup6 results;
			
			 
			if (issparse == false){
				results = solve6<'d', size_t, TFloat>(algorithm, minibatchsize, nthreads, ndata, dimension, data.get(), ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, C_init.get(), data_indices_init_from.get(), setseed, static_cast<size_t>(seed),  static_cast<TFloat>(1000*maxtime), maxiter, voutfn, nvaldata, valdata.get(), valperiod, cmsewritefn, cmserate, gbphi); //24
			}
			else{
				results = solve6<'s', size_t, TFloat>(algorithm, minibatchsize, nthreads, ndata, dimension, &s_data, ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, C_init.get(), data_indices_init_from.get(), setseed, static_cast<size_t>(seed),  static_cast<TFloat>(1000*maxtime), maxiter, voutfn, nvaldata, &s_valdata, valperiod, cmsewritefn, cmserate, gbphi);
			}
	

			//TFloat startmse = std::get<6>(results);
			TFloat endmse = std::get<5>(results);
			size_t niterations = std::get<4>(results);
			size_t duration = std::get<3>(results);
			
			if (moutfn.compare("") != 0){
				std::ofstream mrfile;
				if (ri == 0){
					mrfile.open(moutfn);
					mrfile << "niterations\tduration\tendmse\n";
				}
				else{
					mrfile.open(moutfn, std::ofstream::app);
				}
				mrfile << niterations << "\t" << duration << "\t" << endmse << "\n";
				mrfile.flush();
			}
					
			if (endmse < bestmse){
				/* C, L, inds0, duration, niterations, mse */	
				/* std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<size_t []>, std::unique_ptr<size_t []>, size_t, size_t, TFloat> */
				if (coutfn.compare("") != 0){
					auto C_final = std::move(std::get<0>(results));
					file.open(coutfn);
					size_t index = 0;
					while (index < ncentroids*dimension){
						file << C_final[index] << "\t";
						if ((index != 0) && (index % dimension  == 0)){
							file << "\n";
						}
						++index;
					}
					file.close();
				}
				
	
				if (loutfn.compare("") != 0){
					auto L_final = std::move(std::get<1>(results));
					file.open(loutfn);
					size_t index = 0;
					while (index < ndata){
						file << L_final[index] << "\n";
						++index;
					}
					file.close();
				}
				
				if (ioutfn.compare("") != 0){
					auto i0_final = std::move(std::get<2>(results));
						
					//check that i0_final is not nullptr
					if (i0_final.get() == nullptr){
						file << "required pointer is a null pointer, are you sure that it makes sense to print the initialisation indices to file?";
					}
					else{
						file.open(ioutfn);
						size_t index = 0;
						while (index < ncentroids){
							file << i0_final[index] << "\n";
							++index;
						}
					}
					file.close();
				}
			}
		}
	
	}
	
	//}
	

	void solvewrited(const std::string & algorithm, bool issparse, size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir, const std::string & cinfn, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize, std::string & cmsewritefn, size_t cmserate, double gbphi){
		solvewrite<double>(algorithm, issparse, nruns, nthreads, cout_verbosity, file_verbosity, datainfn, coutfn, loutfn, ioutfn, soutfn, voutfn, moutfn, moutdir, cinfn, ind0fn, init0, setseed, seed, ncentroids, maxiter, maxtime, valinfn, valperiod, minibatchsize, cmsewritefn, cmserate, gbphi); //28
	}
	
	void solvewritef(const std::string & algorithm, bool issparse, size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir, const std::string & cinfn, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize, std::string & cmsewritefn, size_t cmserate, float gbphi){
		solvewrite<float>(algorithm, issparse, nruns, nthreads, cout_verbosity, file_verbosity, datainfn, coutfn, loutfn, ioutfn, soutfn, voutfn, moutfn, moutdir, cinfn, ind0fn, init0, setseed, seed,  ncentroids, maxiter, maxtime, valinfn, valperiod, minibatchsize, cmsewritefn, cmserate, gbphi);
	}


	template <typename TFloat>
	void v_solveioless(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const TFloat * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const TFloat * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, TFloat maxtime, size_t maxrounds, TFloat * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, TFloat & mse, size_t minibatchsize, size_t nvaldata, const TFloat * const valdata, size_t valperiod, bool captureverbose, std::string & verbosestring){
		
		auto tup = solveioless<TFloat>(algorithm, nthreads, ndata, dimension, data,  ncentroids, cout_verbosity, initialisation_method, C_init,  data_indices_init_from, setseed, seed, maxtime, maxrounds, minibatchsize, nvaldata, valdata,  valperiod, captureverbose);
		

		std::memcpy(C, std::get<0>(tup).get(), sizeof(TFloat)*ncentroids*dimension);
		std::memcpy(L, std::get<1>(tup).get(), sizeof(size_t)*ndata);
		if (std::get<2>(tup).get()){
			std::memcpy(inds0, std::get<2>(tup).get(), sizeof(size_t)*ncentroids);
		}
		duration = std::get<3>(tup); 
		niterations = std::get<4>(tup);
		
		mse = std::get<5>(tup);
		
		verbosestring = std::move(std::get<6>(tup));
	
	}
	
	void v_solveiolessf(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const float * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const float * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, float maxtime, size_t maxrounds, float * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, float & mse, size_t minibatchsize, size_t nvaldata, const float * const valdata, size_t valperiod, bool captureverbose, std::string & verbosestring){
		v_solveioless<float>(algorithm, nthreads, ndata, dimension, data, ncentroids, cout_verbosity,initialisation_method,  C_init, data_indices_init_from,  setseed,  seed,  maxtime,  maxrounds, C,  L, inds0, duration, niterations, mse, minibatchsize, nvaldata, valdata, valperiod, captureverbose, verbosestring);		
	}
	
	void v_solveiolessd(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const double * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const double * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, double maxtime, size_t maxrounds, double * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, double & mse, size_t minibatchsize, size_t nvaldata, const double * const valdata, size_t valperiod, bool captureverbose,   std::string & verbosestring){
		v_solveioless<double>(algorithm, nthreads, ndata, dimension, data, ncentroids, cout_verbosity,initialisation_method,  C_init, data_indices_init_from,  setseed,  seed,  maxtime,  maxrounds, C,  L, inds0, duration, niterations, mse, minibatchsize, nvaldata, valdata, valperiod, captureverbose, verbosestring);					
	}
	

}


