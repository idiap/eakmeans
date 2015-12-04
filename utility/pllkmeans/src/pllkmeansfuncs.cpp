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

#include "pllcluster.h"
#include "pllkmeansfuncs_nonvoid.h"
#include "pllkmeansfuncs_void.h"
#include "pllkmeansfuncs.hpp"
#include <iostream>
#include <exception>
	
namespace cluster{
	
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
	void solvewrite(const std::string & algorithm, size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir, const std::string & cinf, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize){
		
		size_t nvaldata = 0;
		std::unique_ptr<TFloat []> valdata;
		
		/* load data */
		size_t ndata;
		size_t dimension;		
		std::fstream dfile(datainfn, std::ios_base::in);
		dfile >> ndata;
		dfile >> dimension;
		std::unique_ptr<TFloat []> data (new TFloat [ndata*dimension]);
		TFloat * raw_dataptr = data.get();
		while (dfile >> *raw_dataptr){
			++raw_dataptr;
    }
		dfile.close();
		
		
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
			std::fstream init_file(ind0fn, std::ios_base::in);
			size_t * rp = data_indices_init_from.get();
			while (init_file >> *rp){
				++rp;
			}
			init_file.close();	
		}
		
		else {
			
		}
		
		if (valinfn.compare("") != 0){
			size_t dimvaldata;		
			std::fstream valdfile(valinfn, std::ios_base::in);
			valdfile >> nvaldata;
			valdfile >> dimvaldata;
			if (dimvaldata != dimension){
				throw std::runtime_error("The dimension of the training data and validation data do not seem to have the same dimension. Recall the format of the data files : first line is (ndata dimension) followed by ndata lines with dimension floating point values"); 
			}
			
			valdata.reset(new TFloat [nvaldata*dimension]);
			TFloat * raw_valdataptr = valdata.get();
			while (valdfile >> *raw_valdataptr){
				++raw_valdataptr;
			}
			valdfile.close();	
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
			
			
			std::cout << "entering pllclustering function" << std::endl;
			std::cout << nvaldata << " " << valdata.get() << " " << valperiod << std::endl;
			
			//std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat> results;
			
			//if (algorithm.compare("minibatch") != 0){
			auto results = solve6<size_t, TFloat>(algorithm, minibatchsize, nthreads, ndata, dimension, data.get(), ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, C_init.get(), data_indices_init_from.get(), setseed, static_cast<size_t>(seed),  static_cast<TFloat>(1000*maxtime), maxiter, voutfn, nvaldata, valdata.get(), valperiod);
			//}
			
			//else{
				//results = solve6<size_t, TFloat>(algorithm, minibatchsize, nthreads, ndata, dimension, data.get(), ncentroids, cout_verbosity, file_verbosity, file, initialisation_method, C_init.get(), data_indices_init_from.get(), setseed, static_cast<size_t>(seed),  static_cast<TFloat>(1000*maxtime), maxiter, voutfn, nvaldata, valdata.get(), valperiod);
			//}
			
		
			//TFloat startmse = std::get<6>(results);
			TFloat endmse = std::get<5>(results);
			size_t niterations = std::get<4>(results);
			size_t duration = std::get<3>(results);
			
			//std::cout << "Final mse : " << runmse << std::endl;
			
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
		
		
			//if (file.is_open()){
				//file << "\nfinal mse: " << runmse  << "\n";
				//file.close();
			//}
			
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
	

	void solvewrited(const std::string & algorithm, size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir, const std::string & cinfn, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize){
		solvewrite<double>(algorithm, nruns, nthreads, cout_verbosity, file_verbosity, datainfn, coutfn, loutfn, ioutfn, soutfn, voutfn, moutfn, moutdir, cinfn, ind0fn, init0, setseed, seed, ncentroids, maxiter, maxtime, valinfn, valperiod, minibatchsize);
	}
	
	void solvewritef(const std::string & algorithm, size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir, const std::string & cinfn, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize){
		solvewrite<float>(algorithm, nruns, nthreads, cout_verbosity, file_verbosity, datainfn, coutfn, loutfn, ioutfn, soutfn, voutfn, moutfn, moutdir, cinfn, ind0fn, init0, setseed, seed,  ncentroids, maxiter, maxtime, valinfn, valperiod, minibatchsize);
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
		
		verbosestring = std::move(std::get<6>(tup));
	
	}
	
	void v_solveiolessf(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const float * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const float * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, float maxtime, size_t maxrounds, float * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, float & mse, size_t minibatchsize, size_t nvaldata, const float * const valdata, size_t valperiod, bool captureverbose, std::string & verbosestring){
		v_solveioless<float>(algorithm, nthreads, ndata, dimension, data, ncentroids, cout_verbosity,initialisation_method,  C_init, data_indices_init_from,  setseed,  seed,  maxtime,  maxrounds, C,  L, inds0, duration, niterations, mse, minibatchsize, nvaldata, valdata, valperiod, captureverbose, verbosestring);		
	}
	
	void v_solveiolessd(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const double * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const double * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, double maxtime, size_t maxrounds, double * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, double & mse, size_t minibatchsize, size_t nvaldata, const double * const valdata, size_t valperiod, bool captureverbose,   std::string & verbosestring){
		v_solveioless<double>(algorithm, nthreads, ndata, dimension, data, ncentroids, cout_verbosity,initialisation_method,  C_init, data_indices_init_from,  setseed,  seed,  maxtime,  maxrounds, C,  L, inds0, duration, niterations, mse, minibatchsize, nvaldata, valdata, valperiod, captureverbose, verbosestring);					
	}
}
