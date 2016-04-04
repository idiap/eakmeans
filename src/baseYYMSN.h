// UNDER CONSTRUCTION.

#ifndef PLL_PLLYINYANGMSNBASEKMEANS_H
#define PLL_PLLYINYANGMSNBASEKMEANS_H

#include "baseYY.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace kmeans{


template <typename TInt, typename TFloat>
/* " max sum norm" for lower bounds versions inherit from here */ 
class YYMSNBase : public kmeans::BaseYY<TInt, TFloat>{
	
	private:
		std::unique_ptr<TFloat []> delta_C;
		std::unique_ptr<TFloat []> u_delta_G;

	protected:	
		TFloat * const get_glowers(){
			return this->get_glowers_base();
		}
		
		TFloat * const get_upb(){
			return this->get_upb_base();
		}
		
		TFloat * const get_delta_C(){
			return delta_C.get();
		}
		
		TFloat * const get_u_delta_G(){
			return u_delta_G.get();
		}
		

		
	public:
		typedef kmeans::BaseYY<TInt, TFloat> YYB;
		template<typename... Args>
		YYMSNBase(Args&&... args): YYB(std::forward<Args>(args)...), 
		delta_C{ new TFloat [this->getncentroids()] },
		delta_G{ new TFloat [this->get_ngroups()] }
		
		{
			this->setalgname("YYMSNBase");
		}
		
		virtual ~YYMSNBase(){}


		virtual TInt get_approximate_memory_requirement(){
			return YYB::get_approximate_memory_requirement() + 
			sizeof(TFloat)*(
				this->getncentroids() + //delta_C
				this->get_ngroups()); //delta_G
		}

		virtual void verbose_write_additional(){
			this->get_verbose_file() << "\n\n ..not implemented down to YYMSNBase..\n\n";
		}

		
		
		
		virtual void set_initialisation_tasks() = 0;

		
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
};

}

#endif

