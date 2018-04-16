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

#ifndef PLLKMEANSFUNCS_HPP
#define PLLKMEANSFUNCS_HPP

#include <sstream>


#include "pllcluster.h"


namespace cluster{
	
	//boilerplate
	template <typename TFloat>
	std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<size_t []>, std::unique_ptr<size_t []>, size_t, size_t, TFloat, std::string>
	
	solveioless(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const TFloat * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const TFloat * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, TFloat maxtime, size_t maxrounds, size_t minibatchsize, size_t nvaldata, const TFloat * const valdata, size_t valperiod, bool captureverbose){
		
    
	
		std::stringstream buffer;		
		auto cout_buff = std::cout.rdbuf();
    
    if (captureverbose == true){
      auto bizzle = buffer.rdbuf();
			std::cout.rdbuf(bizzle);	      
		}	
    
		std::ofstream nowhere;
		
		
		//I assume cmse not wanted
		size_t cmserate = 0;


		TFloat gbphi = 1e-3; //ooph......
    
    
		auto pretro = solve6<'d', size_t, TFloat>(algorithm, minibatchsize, nthreads, ndata, dimension, data, ncentroids, cout_verbosity, 0,  nowhere, initialisation_method, C_init, data_indices_init_from, setseed, seed, maxtime, maxrounds, "", nvaldata, valdata, valperiod, "", cmserate, gbphi);
		

		std::string text;
		
		if (captureverbose == true){
			text = buffer.str();
			std::cout.rdbuf(cout_buff);
		}
		
		else{
			text = "captureverbose was false, so nothing here";
		}
		

		auto retro = std::move(std::tuple_cat(std::move(pretro), std::make_tuple(text)));//, std::make_tuple<std::string>("bwerlk"));		

		return retro;
	//}	
	}	
}

#endif
