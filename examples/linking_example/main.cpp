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

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>

#include "pllkmeansfuncs_nonvoid.h"



int main(){
  
  
  /* full list of algorithms can be found in pllcluster.h, for exact k-means on dense data, you'll probably want one of:
   * exp-ns (exponion - ns), selk-ns (simplified elkan - ns), syin-ns (simplified yinyang - ns) 
   * See ICML paper Fast K-Means with Accurate Bounds, or --help flag in executable or python function string for more info.
   * */
  std::string algorithm = "exp-ns";
  
  size_t nthreads = 2;
  
  /* we make some data */
  size_t ndata = 10000;
  size_t dimension = 3;
  std::vector<float> v_data (ndata*dimension, 0);
  for (size_t i = 0; i < ndata*dimension; ++i){
    v_data[i]= static_cast<float> (rand() % 100);
  }
  
  size_t ncentroids = 100;
  int cout_verbosity = 2;
  std::string initialisation_method = "from_indices";
  const float * const C_init = nullptr;
  std::vector<size_t> v_data_indices_init_from(ncentroids);
  std::iota(v_data_indices_init_from.begin(), v_data_indices_init_from.end(), 0);
  bool setseed = true;
  size_t seed = 1011;
  float maxtime = 1000;
  size_t maxrounds = 1000;
  size_t minibatchsize = 0;
  bool capture_verbose = false;
  
   
  std::cout << "entering solveiolessf" << std::endl;  
  /* the double version is sloveiolessd : see pllkmeansfuncs_nonvoid */
  auto results = cluster::solveiolessf(algorithm, nthreads, ndata, dimension, v_data.data(), ncentroids, cout_verbosity, initialisation_method, C_init,  v_data_indices_init_from.data(),  setseed,  seed,  maxtime, maxrounds, minibatchsize,  capture_verbose);
  
  /* return : C, L, inds0, duration, niterations, mse */
  float * C = std::get<0> (results).get();
  size_t * labels = std::get<1> (results).get();
  size_t * starting_indices_returned = std::get<2> (results).get();
  size_t duration = std::get<3> (results);
  size_t niterations = std::get<4> (results);
  double mse = std::get<5> (results);

  std::vector<size_t> counts (ncentroids, 0);
  for (size_t i = 0; i < ndata; ++i){
    ++counts[labels[i]];
  }
  
  std::cout << "- -- -  - -  --- -  -  -- - - -   --- -  -- -   -- -" << std::endl;
  for (size_t k = 0; k < ncentroids; ++k){
    std::cout << "in cluster " << k << " : " << counts[k] << std::endl;
  }
  
  
  
  return 0;
}
