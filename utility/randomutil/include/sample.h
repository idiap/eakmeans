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

#ifndef SAMPLE_H
#define SAMPLE_H

#include <exception>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <random>
#include <algorithm>


namespace randomutil{ 
namespace sample{

/*TODO based on code at: //codegolf.stackexchange.com/questions/4772/random-sampling-without-replacement
 * but changed to make uninclusive of upperbound!
g by universal ref?? I really need to clarify when universal ref should be used..
see https://isocpp.org/blog/2012/11/universal-references-in-c11-scott-meyers */

//O(max - min) algorithm for almost uniform sampling.

template<typename OutputIterator, typename IntegerType, class URNG>


void range_no_replacement(OutputIterator out, IntegerType n, IntegerType min, IntegerType max, URNG && g){
  if (n < 0)
    throw std::runtime_error("negative sample size");
  if (max < min)
    throw std::runtime_error("invalid range");
  if (n > max-min+1)
    throw std::runtime_error("sample size larger than range");

  while (n>0)
  {
    double r = g()/(RAND_MAX+1.0);
    if (r*(max-min) < n)
    {
      *out++ = min;
      --n;
    }
    ++min;
  }
}


template<typename OutputIterator, typename IntegerType>
void range_no_replacement(OutputIterator out, IntegerType n, IntegerType min, IntegerType max){
  //std::minstd_rand0 generator (time(NULL));
	//range_no_replacement(out, n, min, max, generator);
  range_no_replacement(out, n, min, max, rand);
}




template<typename IntegerType>
std::vector<IntegerType> get_range_no_replacement(IntegerType n, IntegerType min, IntegerType max){
  std::vector<IntegerType> samples(n);
  range_no_replacement(samples.data(), n, min, max);
  return samples;
}




template<typename IntegerType, class URNG>
std::vector<IntegerType> get_range_no_replacement(IntegerType n, IntegerType min, IntegerType max, URNG && g){
  std::vector<IntegerType> samples(n);
  range_no_replacement(samples.data(), n, min, max, g);
  return samples;
}



template<typename IntegerType, class URNG>
std::vector<IntegerType> get_permuted_range(IntegerType n, URNG && g){
  std::vector<IntegerType> shuffled(n, 0);
  std::iota(shuffled.begin(), shuffled.end(), 0);
  std::random_shuffle(shuffled.begin(), shuffled.end(), [&g](IntegerType i){return g()%i;});
  return shuffled;
}

}
}



#endif
