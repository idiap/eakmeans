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

#include <string>
#include <vector>
#include <iostream>

#include "stringutilbase.h"
#include "stringutilclustering.h"


int main(){
	std::string tosplit("boygirl12boy3091843girl458boy90358034girl58boy0girl30974boy98girlboy64girls");
	std::string delim("boy");
	
	auto v = stringutil::split(tosplit, delim);
	
	std::cout << "\nto split : " << tosplit << std::endl;
	std::cout << "\ndelim : " << delim << std::endl;
	std::cout << "\npost split : " <<  std::endl;
	
	for (auto &  x : v){
		std::cout << x << std::endl;
	}
	
	return 1;
}
