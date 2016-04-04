#ifndef PLL_BASEELKANKMEANS__H
#define PLL_BASEELKANKMEANS__H

namespace kmeans{

template <typename TInt, typename TFloat>
class BaseElkan : public kmeans::BaseExact<TInt, TFloat>{
		
	protected:
			
		
	public:
		typedef kmeans::BaseExact<TInt, TFloat> BC;
		template<typename... Args>
		BaseElkan(Args&&... args): BC(std::forward<Args>(args)...)
		
		{
			this->assignmemory_elkan_upper_lowers();
			this->setalgname("elkan base");
		}
		
		virtual ~BaseElkan(){}

		virtual void verbose_write_additional() override {}
		virtual void set_initialisation_tasks() = 0;
		virtual void set_C_tasks() = 0;
		virtual void set_X_tasks() = 0;
		
		virtual TInt get_approximate_memory_requirement(){
			return BC::get_approximate_memory_requirement() + this->get_elkan_base_memory();
		}
};

}

#endif



