#ifndef PLL_ELKANKMEANS_3V0_H
#define PLL_ELKANKMEANS_3V0_H

#include "baseelkan.h"
#include "alg_X_selkSN.h"

namespace kmeans{

/* discrepency in ndcalcs as compared to a3v0 due to not computing CC initially (I propose) */

template <typename TInt, typename TFloat>
class P3V0 : public kmeans::BaseElkan<TInt, TFloat>{
				
	protected:
		TFloat * const get_lowers(){
			return this->elkan_lowers_base.get();
		}
		
		TFloat * const get_upbs(){
			return this->elkan_upper_base.get();
		}
		
		TFloat * const get_delta_C(){
			return this->elkan_delta_C.get();
		}
		
		std::function<void(TInt)> update_L_lowers_upper_S_H_3v0_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				this->pll_principal_X(update_L_lowers_upper_S_H_3v0<TInt, TFloat>, ti, this->get_delta_C(), this->get_L() + x0,  this->get_lowers() + x0*this->getncentroids(), this->get_upbs() + x0, this->round);
			};
		}
		
		
	public:
		typedef kmeans::BaseElkan<TInt, TFloat> EB;
		template<typename... Args>
		P3V0(Args&&... args): EB(std::forward<Args>(args)...)


		{
			this->setalgname("p3v0");
			this->elkan_delta_C.reset(new TFloat [this->getncentroids()]);
		}
		virtual ~P3V0(){}

		virtual TInt get_approximate_memory_requirement(){
			return EB::get_approximate_memory_requirement() + 
			sizeof(TFloat)*this->getncentroids(); // delta_C  
		}

		virtual void verbose_write_additional(){
			this->EB_verbose_write_additional();
			/* anything else to add ? */
		}

		virtual void set_initialisation_tasks(){
			/* all Elkan variants have same initialisation tasks */
			this->ElkBase_set_initialisation_tasks();
		}
	
		virtual void set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX)
			};
		}
		
		virtual void set_X_tasks(){
			this->X_tasks = {
				this->update_L_lowers_upper_S_H_3v0_ati()
			};
		}
};

}

#endif
