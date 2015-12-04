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

#ifndef PLLKMEANSVOIDFUNCS_H
#define PLLKMEANSVOIDFUNCS_H

namespace cluster {
	
	
	

	/* As per nonvoid versions, but C, L, inds0, duration, niterations, mse set inplace. They should be initialised to be the right dimension before entering.  (useful function for Cython so that no messing around with smart pointers, although apparently it is straightforward...) */
	void v_solveiolessf(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const float * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const float * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, float maxtime, size_t maxrounds, float * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, float & mse, size_t minibatchsize, size_t nvaldata, const float * const valdata, size_t valperiod  , bool captureverbose, std::string & verbosestring);


	void v_solveiolessd(const std::string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const double * const data, size_t ncentroids, int cout_verbosity, const std::string & initialisation_method, const double * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, double maxtime, size_t maxrounds, double * const C, size_t * const L, size_t * const inds0,size_t & duration, size_t & niterations, double & mse, size_t minibatchsize, size_t nvaldata, const double * const valdata, size_t valperiod  , bool captureverbose, std::string & verbosestring);		


	/* functions used in for kmeans executable */
	void solvewrited(const std::string & algorithm, size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir,  const std::string & cinf, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize);
	
	void solvewritef(const std::string & algorithm,  size_t nruns, size_t nthreads, int cout_verbosity, int file_verbosity, const std::string & datainfn, const std::string & coutfn, const std::string & loutfn,  const std::string & ioutfn, const std::string & soutfn, const std::string & voutfn, const std::string & moutfn, const std::string & moutdir,  const std::string & cinf, const std::string & ind0fn, const std::string & init0, bool setseed, size_t seed, size_t ncentroids, size_t maxiter, double maxtime, const std::string & valinfn, size_t valperiod, size_t minibatchsize);
	
	
	

}

#endif
