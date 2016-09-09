/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_PLLCLUSTER_H
#define PLL_PLLCLUSTER_H

#include <string>
#include <new>
#include <iostream>
#include <fstream>
#include <cstdlib>


#include "arrutilv2l0.h"

//algorithms with extern template class instantiations in algorithm header (at release time, this should be empty?)
//...

//////algorithms which will be built from scratch here
#include "sparsesimple.h"
#include "sparseelkan3v0.h"
#include "sparseminibatch.h"
#include "sparsestandardminibatch.h"
#include "SparseGBMseSimple.h"
#include "SparseGBMse3v1.h"
#include "GBPSimple.h"
#include "GBMseSimple.h"
#include "GBMse3v1.h"
#include "minibatch.h"
#include "standardminibatch.h"
#include "mb3v0.h"
#include "exactsimplebatch.h"
#include "simple1.h"
#include "elkan3v0.h"
#include "elkan4v2.h"
#include "elkan5v1.h"
#include "elkan6v0.h"
#include "hamerly11v0.h"
#include "hamerly12v6.h"
#include "hamerly12v7.h"
#include "hamerly13v0.h"
#include "YY17v6.h"
#include "YY17v5.h"
#include "YY17v3.h"
#include "YY17v2.h"
#include "YY21v3.h"
#include "YY21v4.h"
#include "YY21v5.h"

//dangerous for some or other reason
//#include "simplest.h"



namespace cluster{
	
	
	typedef unsigned short tautype;
		
	
	template <typename TInt, typename TFloat, typename... Args>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat>
	get_tup6_sparse(const std::string & ENTERING_OPENBLAS_NUM_THREADS, const std::string & algorithm, TInt minibatchsize, Args&&... args)
	{
		
		
		arrutilv2::proxy_openblas_set_num_threads(1);
		std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat> tup;

	

		if (algorithm.compare("sparsesimple") == 0){
			kmeans::SparseSimple<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("sparseselkSN") == 0 || algorithm.compare("sparsep3v0") == 0){
			kmeans::SP3V0<TInt, TFloat> akmeans(std::forward<Args>(args)...);			
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("sparseminibatch") == 0){
			kmeans::SparseMiniBatch<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("sparsestandardminibatch") == 0){
			kmeans::SparseStandardMiniBatch<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}

		else if (algorithm.compare("sparsegbmsesimple") == 0){
			kmeans::SparseGBMseSimple<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}

		else if (algorithm.compare("sparsegbmse3v1") == 0){
			kmeans::SparseGBMse3v1<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}				
		
		else{

			std::string error_string = std::string("Unrecognised algorithm in pllcluster : ") + algorithm + ". ";
			std::cerr << "error_string" << "\nwill now attempt to clean up openblas threads, using ENTERING_OPENBLAS_NUM_THREADS (" << ENTERING_OPENBLAS_NUM_THREADS << ")" << std::endl; 
			if (ENTERING_OPENBLAS_NUM_THREADS.compare("") != 0){
				arrutilv2::proxy_openblas_set_num_threads(stoi(ENTERING_OPENBLAS_NUM_THREADS));
			}
			throw std::runtime_error(error_string);
		}
		
		setenv("OPENBLAS_NUM_THREADS", ENTERING_OPENBLAS_NUM_THREADS.c_str(),1); 
		return tup;
	}
	
	
	template <typename TInt, typename TFloat, typename... Args>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat>
	get_tup6_dense(const std::string & ENTERING_OPENBLAS_NUM_THREADS, const std::string & algorithm, TInt minibatchsize, Args&&... args)
	{
		
		arrutilv2::proxy_openblas_set_num_threads(1);
		std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat> tup;
	
  
	      
		//TODO : change name : minibatchsize 
		if (algorithm.compare("gbsimple") == 0){
			kmeans::GBPSimple<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...); //1, 21
			return akmeans.get_6();
		}
		
		else if (algorithm.compare("gbmse") == 0){
			kmeans::GBMseSimple<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			return akmeans.get_6();
		}
		
		else if (algorithm.compare("gbmse3v1") == 0){
			kmeans::GBMse3v1<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			return akmeans.get_6();
		}
		
		else if (algorithm.compare("minibatch") == 0){
			kmeans::MiniBatch<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		
		else if (algorithm.compare("standardminibatch") == 0){
			kmeans::StandardMiniBatch<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("mb3v0") == 0 || algorithm.compare("elkminibatch") == 0){
			kmeans::MB3V0<TInt, TFloat> akmeans(minibatchsize, std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("exactsimplebatch") == 0 || algorithm.compare("stab") == 0){
			kmeans::SimpleExactBatchKmeans<TInt, TFloat> sokmo(std::forward<Args>(args)...);
			tup = sokmo.get_6();
		}
		
		else if (algorithm.compare("p4v2") == 0 || algorithm.compare("selkNS") == 0  || algorithm.compare("selk-ns") == 0){
			kmeans::P4V2<TInt, tautype, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}

		else if (algorithm.compare("p12v7") == 0 || algorithm.compare("expNS") == 0 || algorithm.compare("exp-ns") == 0){
			kmeans::P12V7<TInt, tautype, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("p21v3") == 0 || algorithm.compare("syinNS") == 0 || algorithm.compare("syin-ns") == 0){
			kmeans::P21V3<TInt, tautype, TFloat> akmeans(std::forward<Args>(args)...);
			return akmeans.get_6();
		}			


		//else if (algorithm.compare("simplest") == 0){
			//kmeans::SimplestKmeans<TInt, TFloat> sokmo(std::forward<Args>(args)...);
			//tup = sokmo.get_6();
		//}
						
		else if (algorithm.compare("simple") == 0 || algorithm.compare("sta") == 0){
			kmeans::SimpleKmeans1<TInt, TFloat> sokmo(std::forward<Args>(args)...);
			tup = sokmo.get_6();
		}
		
		
		else if (algorithm.compare("p3v0") == 0 || algorithm.compare("selkSN") == 0 || algorithm.compare("selk-sn") == 0){
			kmeans::P3V0<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		
		else if (algorithm.compare("p5v1") == 0  || algorithm.compare("elkSN") == 0 || algorithm.compare("elk-sn") == 0){
			kmeans::P5V1<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("p6v0") == 0 || algorithm.compare("elkNS") == 0 || algorithm.compare("elk-ns") == 0){
			kmeans::P6V0<TInt, tautype, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}
		
		else if (algorithm.compare("p11v0") == 0 || algorithm.compare("ham") == 0){
			kmeans::P11V0<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}	

		else if (algorithm.compare("p12v6") == 0 || algorithm.compare("expSN") == 0 || algorithm.compare("exp-sn") == 0){
			kmeans::P12V6<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}	
		
		
		else if (algorithm.compare("p13v0") == 0 || algorithm.compare("ann") == 0 ){
			kmeans::P13V0<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup = akmeans.get_6();
		}	
		
			
		else if (algorithm.compare("p17v6") == 0 || algorithm.compare("yin") == 0){
			kmeans::P17V6<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup= akmeans.get_6();
		}			
				
		else if (algorithm.compare("p17v5") == 0){
			kmeans::P17V5<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			tup =  akmeans.get_6();
		}			
		
				
		else if (algorithm.compare("p17v3") == 0 || algorithm.compare("syinSN") == 0){
			kmeans::P17V3<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			return akmeans.get_6();
		}				
		
		else if (algorithm.compare("p17v2") == 0){
			kmeans::P17V2<TInt, TFloat> akmeans(std::forward<Args>(args)...);
			return akmeans.get_6();
		}				
		

		else if (algorithm.compare("p21v4") == 0){
			kmeans::P21V4<TInt, tautype, TFloat> akmeans(std::forward<Args>(args)...);
			return akmeans.get_6();
		}						
		
		else if (algorithm.compare("p21v5") == 0){
			kmeans::P21V5<TInt, tautype, TFloat> akmeans(std::forward<Args>(args)...);
			return akmeans.get_6();
		}						
		

		else{
      std::cout << "invalid algorithm name" << std::endl;
      //std::cout << "ENTERING_OPENBLAS_NUM_THREADS" << "`" << ENTERING_OPENBLAS_NUM_THREADS << "'" << std::endl;
			
      if (ENTERING_OPENBLAS_NUM_THREADS.compare("") != 0){
        arrutilv2::proxy_openblas_set_num_threads(std::stoi(ENTERING_OPENBLAS_NUM_THREADS));
      }
			std::string error_string = std::string("Unrecognised algorithm in pllcluster : ") + algorithm + ". ";
			throw std::runtime_error(error_string);
		}
		
		setenv("OPENBLAS_NUM_THREADS", ENTERING_OPENBLAS_NUM_THREADS.c_str(),1); 
		return tup;

	}		
	
	
	template <char sparsity, typename TInt, typename TFloat>
	class tup6_getter{
		public:
			template <typename... Args>
			std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat>
			get_tup6(const std::string & ENTERING_OPENBLAS_NUM_THREADS, const std::string & algorithm, TInt minibatchsize, Args&&... args);
	};
	
	template <typename TInt, typename TFloat>
	class tup6_getter<'d', TInt, TFloat>{ //dense version
		public:
			template <typename... Args>
			std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat>
			get_tup6(const std::string & ENTERING_OPENBLAS_NUM_THREADS, const std::string & algorithm, TInt minibatchsize, Args&&... args){
				
				return get_tup6_dense <TInt, TFloat> (ENTERING_OPENBLAS_NUM_THREADS, algorithm, minibatchsize, std::forward<Args>(args)...);
			}
	};
	
	template <typename TInt, typename TFloat>
	class tup6_getter<'s', TInt, TFloat>{ //sparse version
		public:
			template <typename... Args>
			std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat>
			get_tup6(const std::string & ENTERING_OPENBLAS_NUM_THREADS, const std::string & algorithm, TInt minibatchsize, Args&&... args){
				return get_tup6_sparse <TInt, TFloat> (ENTERING_OPENBLAS_NUM_THREADS, algorithm, minibatchsize, std::forward<Args>(args)...);
			}
	};
		
	
	/* C, L, inds0, duration, niterations, mse */
	template <char sparsity, typename TInt, typename TFloat, typename...  Args>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat>
	solve6(const std::string & algorithm, TInt minibatchsize, Args&&... args){//1,1,21
		
		std::string ENTERING_OPENBLAS_NUM_THREADS = "";
		if (std::getenv("OPENBLAS_NUM_THREADS")){
			ENTERING_OPENBLAS_NUM_THREADS = std::getenv("OPENBLAS_NUM_THREADS");
		}
    
    
    //std::cout << "about to enter a try" << std::endl;		
		try{
			tup6_getter<sparsity, TInt, TFloat> a_getter;
			return a_getter.get_tup6(ENTERING_OPENBLAS_NUM_THREADS, algorithm, minibatchsize, std::forward<Args>(args)...);//1,1,1,21
		
		}
		
		
		catch (std::bad_alloc& ba){
			
			setenv("OPENBLAS_NUM_THREADS", ENTERING_OPENBLAS_NUM_THREADS.c_str(),1);
			std::cerr << "bad_alloc caught: " << ba.what() << '\n';
			std::cerr << "the attempt to cluster failed. Probably due to a memory allocation failure : the memory demand was too large (ndata/ncentroids/ngroups/n... too large). For now, will return tuple of nullptrs etc\n" << std::endl; 
		}
		
		return 	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []>, TInt, TInt, TFloat > {nullptr, nullptr, nullptr, 0, 0, 0};
	}
	
	/* time, niters, energy */
	template <char sparsity, typename TInt, typename TFloat, typename...  Args>
	std::tuple<TInt, TInt, TFloat >
	getfinalstats(Args&&... args){
		auto tup6 = solve6<sparsity, TInt, TFloat>(std::forward<Args>(args)...);
		return std::make_tuple(std::get<3>(tup6), std::get<4>(tup6), std::get<5>(tup6));
	}


	/* C, L inds0 */
	template <char sparsity, typename TInt, typename TFloat, typename...  Args>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> >
	solve(Args&&... args){
		auto tup6 = solve6<sparsity, TInt, TFloat>(std::forward<Args>(args)...);
		return std::make_tuple(std::move(std::get<0>(tup6)), std::move(std::get<1>(tup6)), std::move(std::get<2>(tup6)));
	}
	




	/* load data (TODO : include test for file validity) */
	template <char sparsity, typename TInt, typename TFloat, typename...  Args>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> , TInt, TInt, TFloat>
	solve6_fromfile(const std::string & datafilename, const std::string & algorithm, TInt nthreads, Args&&... args){
		TInt ndata;
		TInt dimension;		
		std::fstream dfile(datafilename, std::ios_base::in);
		dfile >> ndata;
		dfile >> dimension;
		std::unique_ptr<TFloat []> data (new TFloat [ndata*dimension]);
		TFloat * raw_dataptr = data.get();
		while (dfile >> *raw_dataptr){
			++raw_dataptr;
    }
		dfile.close();
		return solve6<sparsity, TInt, TFloat>(algorithm, nthreads, ndata, dimension, data.get(), std::forward<Args>(args)...);
	}	
		

	template <char sparsity, typename TInt, typename TFloat, typename...  Args>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, std::unique_ptr<TInt []> >
	solve3_fromfile(Args&&... args){
		auto tup6 = solve6_fromfile<sparsity, TInt, TFloat>(std::forward<Args>(args)...);
		return std::make_tuple(std::move(std::get<0>(tup6)), std::move(std::get<1>(tup6)), std::move(std::get<2>(tup6)));
	}
	
	template <char sparsity, typename TInt, typename TFloat, typename...  Args>
	std::tuple<TInt, TInt, TFloat >
	getfinalstats_fromfile(Args&&... args){
		auto tup6 = solve6_fromfile<sparsity, TInt, TFloat>(std::forward<Args>(args)...);
		return std::make_tuple(std::get<3>(tup6), std::get<4>(tup6), std::get<5>(tup6));
	}

}






#endif

