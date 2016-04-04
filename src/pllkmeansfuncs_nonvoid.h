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
