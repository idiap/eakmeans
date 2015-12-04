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

#include "sample.h"
#include "randomarray.h"

#include <memory>

#include <random>
#include <iostream>

int previous_tests(){
	
	
	srand(time(NULL));
	
	unsigned n = 10;
	unsigned min  = 100;
	unsigned max = 200;
	
	//populate an array with samples from range:
	unsigned * const indices  = new unsigned [n];
	
	std::minstd_rand0 generator (time(NULL) + 1);
	randomutil::sample::range_no_replacement(indices, n, min, max, generator);
	for (unsigned i = 0; i < n; ++ i){
		std::cout << indices[i] << " ";
	}
	std::cout << std::endl;
	
	
	
	//using default number generator
	
	for (unsigned h = 0; h < 5; ++h){
		std::cout << "-+-+-+-+-+-+-" << h << std::endl;
		randomutil::sample::range_no_replacement(indices, n, min, max);
	
		for (unsigned i = 0; i < n; ++ i){
			std::cout << indices[i] << " ";
		}
		std::cout << std::endl;
	}
	
	delete [] indices;
	
	
	//get a vector with samples from range:
	std::vector<unsigned> samples = randomutil::sample::get_range_no_replacement(n, min, max);
	for (unsigned i = 0; i < n; ++ i){
		std::cout << samples[i] << " ";
	}
	std::cout << std::endl;
	
	
	//get a vector with samples from range and generator:
	std::cout << "vector with samples from range and generator" << std::endl;
	samples = randomutil::sample::get_range_no_replacement(n, min, max, generator);
	for (unsigned i = 0; i < n; ++ i){
		std::cout << samples[i] << " ";
	}
	std::cout << std::endl;
	
	
	
	
	
	//fill arrays with random values:
	size_t size = 5;
	
	unsigned *  ua = new unsigned [size];
	int * ia = new int [size];
	float * fa = new float [size];
	double * da = new double [size];
	
	randomutil::randomarray::filluniform(size, ua, 0, 10);
	randomutil::randomarray::filluniform(size, ia, -1, 2);
	randomutil::randomarray::filluniform(size, fa, 3, 3.5);
	randomutil::randomarray::filluniform(size, da, 3, 3.5);
	for (size_t i = 0; i  < size ; ++ i){
		std::cout << ua[i] << " ";
	}
	std::cout << std::endl;
	
	for (size_t i = 0; i  < size ; ++ i){
		std::cout << ia[i] << " ";
	}
	std::cout << std::endl;


	for (size_t i = 0; i  < size ; ++ i){
		std::cout << fa[i] << " ";
	}
	std::cout << std::endl;


	for (size_t i = 0; i  < size ; ++ i){
		std::cout << da[i] << " ";
	}
	std::cout << std::endl;

	
	delete [] ua;
	delete [] ia;
	delete [] fa;
	delete [] da;

	return 0;

}


	
int unit_tests(){
	srand(time(NULL));
	typedef unsigned TInt ;
	typedef float TNumber ;
	std::vector<TNumber> options {0,1,2,77.777};
	TInt size_tofill = 20;
	std::cout << options.size();
	std::vector<TNumber> tofill(size_tofill);
	randomutil::randomarray::filluniform(size_tofill, tofill.data(), options);
	for (TInt i = 0; i < size_tofill; ++i){
		std::cout << tofill[i] << " ";
	}
	std::cout<<std::endl;
	
	
	auto sampled = randomutil::randomarray::getuniform<TNumber>(size_tofill, options);
	for (TInt i = 0; i < size_tofill; ++i){
		std::cout << sampled[i] << " ";
	}
	std::cout<<std::endl;

	return 0;
}

int main(){
	unit_tests();
	return 0;
}
