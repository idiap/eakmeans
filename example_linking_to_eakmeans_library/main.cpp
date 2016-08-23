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
  auto results = cluster::solveiolessf(algorithm, nthreads, ndata, dimension, v_data.data(), ncentroids, cout_verbosity, initialisation_method, C_init,  v_data_indices_init_from.data(),  setseed,  seed,  maxtime, maxrounds, minibatchsize,  capture_verbose);
  
  return 0;
}
