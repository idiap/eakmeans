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

#include "blastemplates.h"
#include <iostream>

template <typename TFloatType>
void gemmexample(){
	unsigned a = 2;
	unsigned b= 3;
	unsigned c = 4;
	TFloatType * M1 = new TFloatType[a*b];
	for (unsigned i = 0; i < a*b; ++i){
		M1[i] = i;
	}
	TFloatType * M2 = new TFloatType[b*c];
	for (unsigned i = 0; i < b*c; ++i){
		M2[i] = i;
	}
	TFloatType * M3 = new TFloatType[a*c];	
	wblas::gemm<TFloatType>(CblasRowMajor, CblasNoTrans, CblasNoTrans, a, c, b, 1.0, M1, b, M2, c, 0.0, M3, c); 
	for (unsigned i = 0; i < a; ++i){
		for (unsigned j = 0; j < c; ++j){
			std::cout << M3[i*c + j] << " ";
		}
		std::cout << std::endl;
	}
}


int main(){
	float * A = new float[1000];
	auto b= wblas::dot(5, A + 1, 1, A+2, 1);
	std::cout << b << std::endl;
	//wblas::finetest();
	//gemmexample<float>();	
	return 0;
}
