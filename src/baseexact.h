#ifndef PLL_BASEEXACTKMEANSTRUE_H
#define PLL_BASEEXACTKMEANSTRUE_H

#include "basekmeans.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseExact : public kmeans::BaseKmeans<TInt, TFloat> {
	
	private:
		virtual void set_summaries() {
			this->set_summaries_exact();
		}
		
		virtual void set_mse() override {
			this->mse = arrutilv2::getmeanl22at(this->ncentroids, this->dimension, this->get_C(), this->ndata, this->data, this->get_L(), this->get_C_l22s(), this->get_data_l22s());
		}

	protected:
	
		virtual void set_initialisation_tasks() = 0;
		virtual void set_X_tasks() = 0;		
		virtual void set_C_tasks() = 0;

		template <typename Function, typename... Args>
		void pll_principal_X(const Function & X_updater, TInt ti, Args&&... args){
			this->base_pll_principal_X(static_cast<TInt> (0), this->ndata, X_updater, ti, std::forward<Args>(args)...);
		}
		
		std::function<void(TInt)> set_L_ati(){
			return this->set_L_ati(0, this->ndata);
		}
	
	public:
		template<typename... Args>
		BaseExact(Args&&... args): kmeans::BaseKmeans<TInt, TFloat>(std::forward<Args>(args)...){}
		
		virtual ~BaseExact(){}
	
};

}

#endif
