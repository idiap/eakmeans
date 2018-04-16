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

#ifndef ENDOFROUNDSTRING_H
#define ENDOFROUNDSTRING_H

#include <string>



namespace stringutil{
namespace clustering{
namespace helper{
std::string getstars();
}


		
namespace pll{


namespace exact{

std::string getstartsummary_v1(std::string algname, size_t memory_usage, float mse, float val_mse);
std::string getstartsummary_v2(std::string algname, size_t memory_usage, float mse, float val_mse);
std::string getroundsummary_v1(size_t roundchanges);
std::string getroundsummary_v2(size_t round, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, float val_mse);
std::string getfinalsummary_v1(size_t round, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse);
std::string getfinalsummary_v2(size_t round, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse);

}

namespace minibatch{
std::string getstartsummary_v1(std::string algname, size_t memory_usage, float val_mse);
std::string getstartsummary_v2(std::string algname, size_t memory_usage, float val_mse);

std::string getroundsummary_v1(size_t round, size_t nsubrounds, size_t roundchanges);
std::string getroundsummary_v2(size_t round, size_t nsubrounds, size_t subround, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, float val_mse);

std::string getfinalsummary_v1(size_t round, size_t nsubrounds, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse);
std::string getfinalsummary_v2(size_t round, size_t nsubrounds, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse);



}


namespace growbatch{


std::string getstartsummary_v1(std::string algname, size_t memory_usage, float val_mse);
std::string getstartsummary_v2(std::string algname, size_t memory_usage, float val_mse);

std::string getroundsummary_v1(size_t roundchanges, bool didgrow);
std::string getroundsummary_v2(size_t round, size_t nactive, float d_C__by__d_AB, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, float val_mse);

std::string getfinalsummary_v1(size_t round, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse);
std::string getfinalsummary_v2(size_t round, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse);

}




} //end pll




}
}




#endif
