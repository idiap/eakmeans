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

#ifndef PLL_SIMPLEKMEANS_H
#define PLL_SIMPLEKMEANS_H

#include "basesimpleexact.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class SimpleKmeans1 : public kmeans::BaseSimpleExactKmeans<TInt, TFloat>{


	protected: 
	
			virtual void set_X_tasks(){
		
			this->X_tasks = {

				arrutilv2::update_L_S_H_ati(this->nthreads, this->ndata, this->dimension, this->data, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dsums(), this->get_dcounts(), this->get_sums(), this->get_counts(), this->nchanges, this->work_mutex, this->ndcalcs_X)
			
			};
			
		}

	public:		
		template<typename... Args>
		SimpleKmeans1(Args&&... args): kmeans::BaseSimpleExactKmeans<TInt, TFloat> (std::forward<Args>(args)...) {this->setalgname("simple kmeans");
		
		}		
		virtual ~SimpleKmeans1(){};
		
};

}

//extern template class kmeans::SimpleKmeans1<size_t, double>;
//extern template class kmeans::SimpleKmeans1<size_t, float>;

#endif
			
	//protected:
		//virtual void set_X_tasks(){
			
			//std::cout << "in set X taksss ---- " << std::endl;
				
				
				
							//set Ls

			//this->X_tasks = {

				//arrutilv2::update_L_S_H_ati(this->nthreads, this->ndata, this->dimension, this->data, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dsums(), this->get_dcounts(), this->get_sums(), this->get_counts(), this->nchanges, this->work_mutex, this->ndcalcs_X)
			
			//};
			
			//std::cout << "have mojo " << std::endl;
			
			//this->X_tasks.push_back(mojo);
				//[](TInt i){std::cout << "big" << std::endl;}, [](TInt i){std::cout << "bog" << std::endl;} 
				
			//};
			
			//std::cout << "finsished with X taksss" << std::endl;
			
		//}	




			//this->X_tasks = {
				//[this](TInt ti){
					
					//if (ti == 0){
						//this->nchanges = 0;

						//for (TInt i = 0; i < this->getndata(); ++i){
							//TFloat best_distance = std::numeric_limits<TFloat>::max();
							//TInt oldlabel = this->get_L()[i];
							//for (TInt ci = 0; ci < this->getncentroids(); ++ci){
								//TFloat distance2 = 0;
								//for (TInt di = 0; di < this->getdimension(); ++di){
									//TFloat diffy = this->get_C()[ci*this->getdimension() + di] - this->getdata()[i*this->getdimension() + di];
									//distance2 += diffy*diffy;
								//}
								//TFloat distance = std::sqrt(std::max(static_cast<TFloat>(0), distance2));
								//if (distance <= best_distance){
									//best_distance = distance;
									//this->get_L()[i] = ci;
								//}
							//}
							//if (this->get_L()[i] != oldlabel){
								//this->nchanges += 1;
							//}
						//}
					//}
				//}
			//};
