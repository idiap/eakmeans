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

#ifndef PLLKMEANSFUNCS_H
#define PLLKMEANSFUNCS_H


#include <string>
#include <memory>
#include <tuple>
#include <vector>


namespace cluster{


	/* useful for direct use in C++ code */
	/* return : C, L, inds0, duration, niterations, mse */
	std::tuple<std::unique_ptr<float []>, std::unique_ptr<size_t[]>, std::unique_ptr<size_t[]>, size_t, size_t, float, std::string>
	solveiolessf(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const float * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const float * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, float maxtime, size_t maxrounds, size_t minibatchsize, bool captureverbose);
	
	std::tuple<std::unique_ptr<double []>, std::unique_ptr<size_t[]>, std::unique_ptr<size_t[]>, size_t, size_t, double, std::string>
	solveiolessd(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const double * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const double * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, double maxtime, size_t maxrounds, size_t minibatchsize, bool captureverbose);
	
	
}

#endif
