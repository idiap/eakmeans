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

// See http://en.cppreference.com/w/cpp/algorithm/random_shuffle for inspiration :)
template<typename TInt, typename TFloat>
void inplace_shuffle_by_row(TInt nrows, TInt ncols, TFloat * const data){
  
  std::unique_ptr<TFloat []> ptrtemp ( new TFloat [ncols] );
  TFloat * const temp = ptrtemp.get();
  TInt copyindex;
  
  
  for (int i = nrows-1; i > 0; --i) {
    copyindex = rand()%(i+1);
    
    std::memcpy(temp, data + i*ncols, sizeof(TFloat)*ncols);
    std::memcpy(data + i*ncols, data + copyindex*ncols, sizeof(TFloat)*ncols);
    std::memcpy(data + copyindex*ncols, temp, sizeof(TFloat)*ncols);
  }
}



}
}



#endif
