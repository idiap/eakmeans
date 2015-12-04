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

#ifndef PLL_HAMERLYBASEKMEANS_H
#define PLL_HAMERLYBASEKMEANS_H


namespace kmeans{

template <typename TInt, typename TFloat>
class BaseHamerly : public kmeans::BaseExact<TInt, TFloat>{
	
	private:
	
		std::unique_ptr<TFloat []> CC;
		std::unique_ptr<TFloat []> halfminCC;
		std::unique_ptr<TFloat []> delta_C;
		std::unique_ptr<TFloat []> lower_base;
		std::unique_ptr<TFloat []> upper_base;
		
		
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_lower_upper_mati(){
		
			std::vector<std::function<void(TInt)> > tasks;
			
			tasks = this->makeset_C_C_l22s_inds0_mati();	
			tasks.push_back(	
				[this](TInt ti){
					TInt local_ndcalcs = 0;
					TInt x0 = (ti*this->getndata())/this->getnthreads();
					TInt x1 = ((ti+1)*this->getndata())/this->getnthreads();	
					arrutilv2::set_L2_dn(x1 - x0, this->getdimension(), this->getdata() + x0*this->getdimension(), this->getncentroids(), this->get_C(), this->get_data_l22s() + x0, this->get_C_l22s(), this->get_L() + x0, this->get_upper_base() + x0, this->get_lower_base() + x0, local_ndcalcs);
					this->ndcalcs_notX += local_ndcalcs;
				}
			);
			tasks.push_back(
			//set starting mse. TODO : parallelise.	
				[this](TInt ti){	
					if (ti == 0 && this->get_initialisation_method().compare("kmeans++") != 0){ // && this->getfileptr()->is_open()
						this->mse = 0;
						for (TInt i = 0; i < this->getndata(); ++i){
							this->mse += (this->get_upper_base()[i])*(this->get_upper_base()[i]);
						}
						this->mse /= static_cast<TFloat>(this->getndata());
					}
					
				}
			);

			return tasks;	
		}
			
			
	protected:
	
		TFloat * const get_CC(){
			return CC.get();
		}
		
		TFloat * const get_halfminCC(){
			return halfminCC.get();
		}
		
		TFloat * const get_lower_base(){
			return lower_base.get();
		}
		
		TFloat * const get_upper_base(){
			return upper_base.get();
		}
		
		TFloat * const get_delta_C(){
			return delta_C.get();
		}
	
		std::vector<std::function<void(TInt)> > makeset_C_C_l22s_L_inds0_lower_upper_S_H_mati(){
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_lower_upper_mati();

			
			 
			auto init_task_B = arrutilv2::set_S_H_ati(this->nthreads, this->ndata, this->dimension, this->data, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dsums(), this->get_dcounts(), this->get_sums(), this->get_counts(), this->work_mutex);	

			auto initialisation_tasks = std::move(init_tasks_A);
			initialisation_tasks.push_back(std::move(init_task_B));
			return initialisation_tasks;
		}	
			
	public:
		template<typename... Args>
		BaseHamerly(Args&&... args): kmeans::BaseExact<TInt, TFloat>(std::forward<Args>(args)...), 
		
		CC{ new TFloat [this->getncentroids()*this->getncentroids()] },
		halfminCC{ new TFloat [this->getncentroids()] },
		delta_C{ new TFloat [this->getncentroids()] },
		lower_base{ new TFloat [this->getndata()]  },
		upper_base{ new TFloat [this->getndata()] }
		
		{
			this->setalgname("BaseHamerly");
		}
		virtual ~BaseHamerly(){}

		virtual void verbose_write_additional(){
			kmeans::BaseExact<TInt, TFloat>::verbose_write_additional();
			this->get_verbose_file() << "\nlower_base:\n" << lower_base[0] << "\n";			
			this->get_verbose_file() << "\n\nupper_base:\n" << upper_base[0] << "\n";
			/* anything else to print ? */
		}

		virtual void set_initialisation_tasks() = 0;
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;

};

}

#endif



