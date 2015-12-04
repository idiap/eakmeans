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

#ifndef PLL_BASEELKANKMEANS__H
#define PLL_BASEELKANKMEANS__H

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseElkan : public kmeans::BaseExact<TInt, TFloat>{
	
	private:
		std::unique_ptr<TFloat []> lowers_base;
		std::unique_ptr<TFloat []> upper_base;	
		
	protected:
		TFloat * const get_lowers_base(){
			return lowers_base.get();
		}
		
		TFloat * const get_upper_base(){
			return this->upper_base.get();
		}
		
		
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_lowers_upper_mati(){
	
			std::vector<std::function<void(TInt)> > tasks;
				
				tasks = this->makeset_C_C_l22s_inds0_mati();				
				tasks.push_back(	
				
				[this](TInt ti){
					TInt local_ndcalcs = 0;
					TInt x0 = (ti*this->getndata())/this->getnthreads();
					TInt x1 = ((ti+1)*this->getndata())/this->getnthreads();
					
					arrutilv2::set_rrl2ss_argminmins(x1 - x0, this->getdimension(), this->getdata() + x0*this->getdimension(), this->getncentroids(), this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_lowers_base() + x0*this->getncentroids(), this->get_L() + x0, this->get_upper_base() + x0, local_ndcalcs);
					this->ndcalcs_notX += local_ndcalcs;
					}
					
				);
				
				
				//set mse. TODO : paralellise
				tasks.push_back(
				
				[this](TInt ti){	
					if (ti == 0 && this->get_initialisation_method().compare("kmeans++") != 0){
						this->mse = 0;
						for (TInt i = 0; i < this->getndata(); ++i){
							this->mse += (this->get_upper_base()[i])*(this->get_upper_base()[i]);
						}
						this->mse /= static_cast<TFloat>(this->getndata());
					}
					
				}
				
				);

			//}
			return tasks;
		}
		
		
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati(){

			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_lowers_upper_mati();

			auto init_task_B = this->set_S_H_ati();
	
			auto initialisation_tasks = std::move(init_tasks_A);
			//initialisation_tasks.push_back(std::move(init_task_Af));
			initialisation_tasks.push_back(std::move(init_task_B));
	
			return initialisation_tasks;
		}
		
	
		void EB_set_initialisation_tasks(){
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_lowers_upper_S_H_mati();
		}



		
		
	public:
		typedef kmeans::BaseExact<TInt, TFloat> BC;
		template<typename... Args>
		BaseElkan(Args&&... args): BC(std::forward<Args>(args)...), 
		
		lowers_base{ new TFloat [this->getndata()*this->getncentroids()]  },
		upper_base{ new TFloat [this->getndata()] }
		{
			this->setalgname("elkan base");
		}
		virtual ~BaseElkan(){}

		virtual void verbose_write_additional(){}
		
		void EB_verbose_write_additional(){
		
			this->get_verbose_file() << "\nlowers base:\n";
			for (TInt ci = 0; ci < this->getncentroids(); ++ci){
				this->get_verbose_file() << lowers_base[ci] << "\t";
			}
			this->get_verbose_file() << "\n\nupper base:\n" << upper_base[0] << "\n";
		}

		virtual void set_initialisation_tasks(){}

		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
		
		virtual TInt get_approximate_memory_requirement(){
			return BC::get_approximate_memory_requirement() + 
			sizeof(TFloat)*this->getndata()*this->getncentroids() + //lower base
			sizeof(TFloat)*this->getndata(); //upper base
		}
};

}

#endif



