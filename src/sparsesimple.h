#ifndef PLL_SPARSESIMPLE_H
#define PLL_SPARSESIMPLE_H

#include "basesparseexact.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class SparseSimple : public kmeans::BaseSparseExact<TInt, TFloat> {
	
	private: 
		
	public:
		template<typename... Args>
		SparseSimple(Args&&... args): kmeans::BaseSparseExact<TInt, TFloat> (std::forward<Args>(args)...) {
			this->setalgname("sparse-simple-kmeans");
		}
		
		virtual ~SparseSimple(){};
	
	protected:
		virtual void set_initialisation_tasks(){
			auto init_tasks_A = this->makeset_C_C_l22s_L_inds0_mati();
			auto init_task_B = 
				[this](TInt ti){
					if (ti == 0){
						this->set_S_H(static_cast<TInt>(0), this->ndata);
					}
				};
			this->initialisation_tasks = std::move(init_tasks_A);
			this->initialisation_tasks.push_back(std::move(init_task_B));
		}
		
			
		virtual void set_X_tasks() override final{
			
			this->X_tasks = {
				
				
				//experiments show that this is were the majority of the time is spent
				sparse::update_L_label_changes_ati(this->nthreads, *this->ptrdata, this->ncentroids, this->get_C(), this->get_data_l22s(), this->get_C_l22s(), this->get_L(), this->ndcalcs_X, this->where_label_changes),

			
				[this](TInt ti){
					if (ti == 0){
						sparse::update_S_H_from_label_changes(*this->ptrdata, this->where_label_changes, this->get_sums(), this->get_counts(), this->nchanges);
					}
				},
				
				[this](TInt ti){
					if (ti == 0){
					}
				} 			
				
			};
		}
		
		virtual void set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_from_SH_ati(this->nthreads, this->ncentroids, this->dimension, this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s())
			};
		}
		
		
		virtual void verbose_write_additional(){
				//TODO
		}

};
} 


#endif
