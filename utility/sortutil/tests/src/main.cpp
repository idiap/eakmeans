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

#include <iostream>
#include <vector>
#include <chrono>

#include "sortutil.h"



void argsort_test(){
	std::vector<float>  A {1.1, 4.4, 3.3, 5.5, 0.0};
	auto asorted = sort::get_argsort_decreasing<unsigned, float>(A);		
	auto dsorted = sort::get_argsort_increasing<unsigned, float>(A);		
	for (unsigned i = 0; i < 5; ++i){
		std::cout << asorted[i] << "\t" << A[asorted[i]] << std::endl;
	}
	
	for (unsigned i = 0; i < 5; ++i){
		std::cout << dsorted[i] << "\t" << A[dsorted[i]] << std::endl;
	}
	
	
}

void geometric_elements_test(){
	std::vector<double> data {
		1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
		17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,
		17,14,6,4,3,15,9,8,7,2,1,13,5,12,11,10,16,
		15,9,8,7,10,16,2,17,14,6,4,3,1,13,5,11,12
	};
		
	unsigned nrows = 3;
	unsigned ncols = 17;
	unsigned npartitions = 5;
	sort::geometric_elements(nrows, ncols, npartitions, data.data(), std::less <double> ());	
	
	std::cout << "\n";
	for (unsigned i = 0; i < nrows; ++i){
		for (unsigned c = 0; c < ncols; ++ c){
			std::cout << data[i*ncols + c] << "\t";
		}
		std::cout << std::endl;
	}
}

void speedtest_nthelement(){
	//if almost sorted, much faster. use this fact!
	
	int nvals = 10000000;
	std::vector<int> ordered (nvals);
	std::vector<int> almostordered (nvals);
	std::vector<int> random (nvals);
	
	for (int i = 0; i < nvals; ++i){
		ordered[i] = i;
		almostordered[i] = i + rand()%10;
		random[i] = rand()%nvals;
	}
	
	
	auto t0 = std::chrono::high_resolution_clock::now();
	std::nth_element(ordered.data(), ordered.data() + nvals/2, ordered.data()+nvals);	
	auto t1 = std::chrono::high_resolution_clock::now();
	double t_total = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	std::cout << "time to get median of sorted : " << t_total << std::endl;
	
	
	t0 = std::chrono::high_resolution_clock::now();
	std::nth_element(almostordered.data(), almostordered.data() + nvals/2, almostordered.data()+nvals);	
	t1 = std::chrono::high_resolution_clock::now();
	t_total = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	std::cout << "time to get median of almostsorted : " << t_total << std::endl;


	t0 = std::chrono::high_resolution_clock::now();
	std::nth_element(random.data(), random.data() + nvals/2, random.data()+nvals);	
	t1 = std::chrono::high_resolution_clock::now();
	t_total = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	std::cout << "time to get median of random : " << t_total << std::endl;	
}	
	

void test_update_halfordered_partvals_etc1(){
	unsigned ncentroids = 13;
	unsigned nrows = 3;
	unsigned nparts = static_cast<unsigned>(std::floor(std::log2(static_cast<double>(ncentroids - 1))));
	unsigned sihaCC = 2*(std::pow(2, nparts - 1) - 1);
	std::cout << "ncents : " << ncentroids << "\tnrows : " << nrows << "\tnparts : " << nparts << std::endl;
	
	std::vector<int> CC (ncentroids*nrows);
	std::vector<std::pair<int, unsigned>> hops(ncentroids*nrows);
	std::vector<int> partvals((nparts - 1)*nrows);
	std::vector<unsigned> ihaCC(nrows*sihaCC);
	for (unsigned i = 0; i < ncentroids*nrows; ++i){
		CC[i] = 10 + rand() % 90;
		hops[i].first = CC[i];
		hops[i].second = i %ncentroids;
	}
	sort::update_halfordered_partvals_etc1(ncentroids, nrows, nparts,  CC.data(), hops.data(), partvals.data(), ihaCC.data(), std::less<std::pair<int, unsigned>>());
	for (unsigned r = 0; r < nrows; ++r){
		std::cout << "\n\nCC\n";
		for (unsigned i = 0; i < ncentroids; ++i){
			std::cout << CC[r*ncentroids + i] << " ";
		}
		
		std::cout <<"\nhops\n";
		for (unsigned i = 0; i < ncentroids; ++i){
			std::cout << "(" << hops[r*ncentroids + i].first << " " << hops[r*ncentroids + i].second << ") ";
		}
		
		std::cout <<"\npartvals\n";
		for (unsigned i = 0; i < nparts - 1; ++i){
			std::cout << partvals[r*(nparts - 1) + i] << " ";
		}
		
		std::cout <<"\nihaCC\n";
		for (unsigned i = 0; i < sihaCC; ++i){
			std::cout << ihaCC[r*sihaCC + i] << " ";
		}
	}	
	std::cout << "\n\n";
}


int main(){

	//partition_left_right_test();
	//argsort_test();
	//geometric_elements_test();
	//speedtest_nthelement();
	
	test_update_halfordered_partvals_etc1();
	
	

}




