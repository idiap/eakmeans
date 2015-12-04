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

#include "optionsutil.h"
#include <limits>
#include <exception>
#include <iostream>


int main(int argc, char* argv[]){
	
	optionsutil::Options uopts;
	
	
	
	uopts.add("algorithm", "alg", "The algorithm to use, options are sta, expSN, expNS, selkSN, selkNS, elkSN, elkNS, syinNS, syinSN", "s", "expNS");
	std::string algorithm = uopts.options["algorithm"].defval;

	
	uopts.add("nthreads", "nth", "How many threads to use", "i", "1");
	unsigned nthreads = std::stoi(uopts.options["nthreads"].defval);	
	

	uopts.add("fprecision", "fpr", "Floating point precision, options are 32 (float) and 64 (double) ", "i", "64");
	unsigned fprecision = std::stoi(uopts.options["fprecision"].defval);


	
	
	uopts.add("verbosity", "ver", "Verbosity, options are 0 (silent) 1 (acceptable) 2 (verbose) ", "i", "1");	
	int verbosity = std::stoi(uopts.options["verbosity"].defval);
	
	uopts.add("datainfn", "din", "Absolute path of input data, first line of the file should be nrows ncols, subsequent lines should be data, see example", "s", "");	
	std::string datainfn = uopts.options["datainfn"].defval;

	
	uopts.add("coutfn", "cou", "Absolute path where final centroids should be written", "s", "");	
	std::string coutfn = uopts.options["coutfn"].defval;
	
	uopts.add("loutfn", "lou", "Absolute path where final labels should be written", "s", "");	
	std::string loutfn = uopts.options["loutfn"].defval;	



	uopts.add("cinfn", "cin", "Absolute path from where initial centroids should be read, first line of the file should be ``nrows ncols'', subsequent lines should be centroids, see example", "s", "");	
	std::string cinfn = uopts.options["cinfn"].defval;	

	
	
	uopts.add("seed", "see", "Initialisation seed", "i", "time(NULL)");	
	unsigned seed = time(NULL);
	
	uopts.add("ind0fn", "iin", "Absolute path from where indices of data for initialising centroids should be read", "s", "");	
	std::string ind0fn = uopts.options["ind0fn"].defval;
	

	
	uopts.add("init0", "ini", "How to initialise the data, eventually should be one {uniform, kmeanspp, ...}, but currently one of {uniform} ", "s", "");	
	std::string init0 = uopts.options["init0"].defval;


	uopts.add("maxiter", "max", "The maximum number of iterations before terminating ", "i", "infinity");	
	unsigned maxiter = std::numeric_limits<unsigned>::max();
	
	uopts.add("phistop", "phi", "Proportion of changing labels at which to terminate ", "f", "0.0");	
	double phistop = 0.0;
	
	uopts.tail = "Compulsory options are {datainfn}.  Must have at least 1 of {coutfn, loutfn}. Must have exactly 1 of {cinfn, ind0fn, init0}. Optional are {algorithm, nthreads, fprecision, verbosity, seed, maxiter, phistop}. ";
	

 
	
	try {
		uopts.print(std::stoi(argv[1]), std::stoi(argv[2]));
	}
	
	catch (const std::invalid_argument& ia){
		std::cerr << "pass int int as parameters (we're testing the display of options, tab1 tab2.)" << std::endl;
		return 1;
	}
	
	std::cout << algorithm << nthreads << fprecision << verbosity << coutfn << loutfn << cinfn << seed << ind0fn << init0 << maxiter << phistop << std::endl;
	
	return 0; 
	
}
