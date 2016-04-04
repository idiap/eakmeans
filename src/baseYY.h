/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_YINYANGBASEKMEANS_H
#define PLL_YINYANGBASEKMEANS_H

#include "baseexact.h"
#include "simple1.h" //for initial clustering of the centroids


#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{
	

template <typename TInt, typename TFloat>
class BaseYY : public kmeans::BaseExact<TInt, TFloat>{
	
	private:
		TInt ngroups;
		std::unique_ptr<TInt []> groupparts; //ngroups + 1 [0, 7, 13, ... , ncentroids]
		std::unique_ptr<TInt []> groupsizes; //ngroups [0,7, 6, .... 12]
		std::unique_ptr<TInt []> group; //which group does data point belong to
		std::unique_ptr<TFloat []> glowers_base; //ndata x ngroups
		std::unique_ptr<TFloat []> upb_base;

	protected:
		TInt * get_groupparts(){
			return groupparts.get();
		}
		
		TInt * get_groupsizes(){
			return groupsizes.get();
		}
		
		TInt * get_group(){
			return group.get();
		}
	
	
		TFloat * const get_glowers_base(){
			return glowers_base.get();
		}
		
		TFloat * const get_upb_base(){
			return upb_base.get();
		}
		
		
		
		
		std::vector<std::function<void(TInt)>> makeset_C_C_l22s_L_inds0_groups_glowers_upb_base_S_H_mati(){

			auto init_tasks_A = this->exact_makeset_C_C_l22s_inds0_mati();
			

		

			
			/* set C and C_l22s to be clustered, set groupparts and groupsizes */
			auto init_task_B = [this](TInt ti){
				if (ti == 0){
					this->make_centroids_clustered();
				}
			}; 
				
			auto init_task_C = arrutilv2::set_L_group_dn_glowers_ati(this->getnthreads(), this->getndata(), this->getdimension(), this->getdata(), this->getncentroids(), this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_ngroups(), this->get_groupparts(), this->get_groupsizes(), this->get_L(), this->get_group(), this->get_upb_base(), this->get_glowers_base(), this->ndcalcs_notX);
			
			
			//set starting mse. TODO : parallelise.	
			auto init_task_Cb = [this](TInt ti){	
				if (ti == 0 && this->get_initialisation_method().compare("kmeans++") != 0){
					this->mse = 0;
					for (TInt i = 0; i < this->getndata(); ++i){
						this->mse += (this->get_upb_base()[i])*(this->get_upb_base()[i]);
					}
					this->mse /= static_cast<TFloat>(this->getndata());
				}
				
			};
	
			/* setting S and H (could be joined to previous task for 1 fewer barrier in initialisation) */
			auto init_task_D = arrutilv2::set_S_H_ati(this->getnthreads(), this->getndata(), this->getdimension(), this->getdata(), this->getncentroids(), this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->get_dsums(), this->get_dcounts(), this->get_sums(), this->get_counts(), this->work_mutex);			
			
			auto initialisation_tasks = std::move(init_tasks_A);
			initialisation_tasks.push_back(std::move(init_task_B));
			initialisation_tasks.push_back(std::move(init_task_C));
			initialisation_tasks.push_back(std::move(init_task_Cb));
			initialisation_tasks.push_back(std::move(init_task_D));
			
			return initialisation_tasks;
			
		}
		
		
		/* rearrange C and C_l22s, set groupparts and groupsizes */
		void make_centroids_clustered(){
			
			/* cluster centroids starting from first ngroups centroids */
			std::vector<TInt> range(this->get_ngroups());
			std::iota(range.begin(), range.end(), this->get_ngroups());
			
			std::ofstream deffile;
			
			/* perform clustering of initial centroids, with cout_verbosity that of this */
			kmeans::SimpleKmeans1<TInt, TFloat> sokmo(this->getnthreads(), this->getncentroids(),this->getdimension(), this->get_C(), this->get_ngroups(), this->get_cout_verbosity(), deffile, range.data());
			
			range.resize(0);

			/* get rearranged C, C_l22s, groupparts and groupsizes*/
			auto C_etc = sokmo.get_contig_by_cluster_4(static_cast<TInt>(3)); /* no fewer than 3 per cluster */
			
			/*copy rearranged C, C_l22s to C, C_l22s (ideally would do a move, but this deprivatises access to unique pointers, not good) */
			
			std::memcpy(this->get_C(), std::get<0>(C_etc).get(), this->getncentroids()*this->getdimension()*sizeof(TFloat));
			std::memcpy(this->get_C_l22s(), std::get<1>(C_etc).get(), this->getncentroids()*sizeof(TFloat));
			
			/* copy groupparts and groupsizes (no stealing, although could have been) */

			std::memcpy(groupparts.get(), std::get<2>(C_etc).get(), sizeof(TInt)*(ngroups + 1));
			std::memcpy(groupsizes.get(), std::get<3>(C_etc).get(), sizeof(TInt)*ngroups);

		}
	
		void yinyang_initialisation_tasks(){						
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_groups_glowers_upb_base_S_H_mati();		
		}		
		
		
	public:
		typedef kmeans::BaseExact<TInt, TFloat> BC;
		template<typename... Args>
		BaseYY(Args&&... args): BC(std::forward<Args>(args)...), 
		
		/* groupparts and groupsizes not "stolen" from centroid clustering, could have been but lambda function captures become subtle */
		ngroups{ 1 + (this->getncentroids() - 1)/10 },
		groupparts{ new TInt[ngroups + 1] },
		groupsizes{ new TInt[ngroups] },
		group{ new TInt [this->getndata()] },
		glowers_base{ new TFloat [this->getndata()*this->get_ngroups()] },
		upb_base{ new TFloat [this->getndata()] }

		{
			this->setalgname("BaseYY");
		}
		
		virtual ~BaseYY(){}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to YY..\n\n";
		}

		virtual void set_initialisation_tasks() = 0;
		
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
		
		TInt get_ngroups(){
			return ngroups;
		}
		
		
		virtual TInt get_approximate_memory_requirement(){
			return BC::get_approximate_memory_requirement() +
			sizeof(TInt)*(
			this->get_ngroups() + 1 + //groupparts
			this->get_ngroups() + //groupsizes
			this->getndata()) + //group
			sizeof(TFloat)*(this->getndata()*this->get_ngroups() +  //glowers_base
			this->getndata()); //upb_base
		}
	
};

}

#endif

