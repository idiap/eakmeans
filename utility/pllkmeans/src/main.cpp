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

#include "pllkmeansfuncs_nonvoid.h"
#include "pllkmeansfuncs_void.h"
#include <numeric>
#include <stdlib.h> 




int main (){
	
	
	srand(time(NULL));
	typedef float ftype;
	size_t ndata = 1000;
	size_t dimension = 10;
	
	size_t ncentroids = 100;
	std::vector<ftype> data (ndata*dimension, 1);
	std::iota(data.begin(), data.end(), 1);
	
	//auto bot = cluster::solveiolessf("p12v7", 1, ndata, dimension, data.data(), ncentroids, 2, "kmeans++", nullptr, nullptr, false, 1011, 1e44, 9999);
	
	
	
	std::unique_ptr<ftype [] > C (new ftype [dimension*ncentroids] );
	std::unique_ptr<size_t [] > L (new size_t [ndata] );
	std::unique_ptr<size_t [] > inds0 (new size_t [ncentroids] );
	size_t duration;
	size_t niterations;
	ftype mse;
	
	
	std::string emptystring = "";
	
	cluster::v_solveiolessf("exactsimplebatch", 1, ndata, dimension, data.data(), ncentroids, 2, "kmeans++", nullptr, nullptr, false, 1011, 1e44, 9999
	, C.get(), L.get(), inds0.get(), duration, niterations, mse, 
	
	123123123, 
	
	0, nullptr, 0, false, emptystring
	);
		
	return 0;
}
