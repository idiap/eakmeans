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

/*TODO: all variants */


}




}




}
}




#endif
