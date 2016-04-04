#ifndef PLL_SIMPLESTKMEANS_H
#define PLL_SIMPLESTKMEANS_H

#include <limits>

#include "basekmeans.h"

namespace kmeans{

template <typename TInt, typename TFloat>
class SimplestKmeans : public kmeans::BaseExact<TInt, TFloat>{
	
		
	public:
		template<typename... Args>
		/* variadic args ala Eli Bendersky */
		SimplestKmeans(Args&&... args): kmeans::BaseExact<TInt, TFloat> (std::forward<Args>(args)...) {this->setalgname("simplest");}		
		virtual ~SimplestKmeans(){};
			
	protected:
		virtual void set_initialisation_tasks(){
			
			auto init_tasks_A = kmeans::BaseExact<TInt, TFloat>::makeset_C_C_l22s_L_inds0_mati();
			
			auto init_task_B = kmeans::BaseExact<TInt, TFloat>::set_S_H_ati();
			
			this->initialisation_tasks = std::move(init_tasks_A);
			this->initialisation_tasks.push_back(std::move(init_task_B));
			
		}

		virtual void set_X_tasks(){
			//set Ls
			this->X_tasks = {
				[this](TInt ti){
					
					if (ti == 0){
						this->nchanges = 0;

						for (TInt i = 0; i < this->getndata(); ++i){
							TFloat best_distance = std::numeric_limits<TFloat>::max();
							TInt oldlabel = this->get_L()[i];
							for (TInt ci = 0; ci < this->getncentroids(); ++ci){
								TFloat distance2 = 0;
								for (TInt di = 0; di < this->getdimension(); ++di){
									TFloat diffy = this->get_C()[ci*this->getdimension() + di] - this->getdata()[i*this->getdimension() + di];
									distance2 += diffy*diffy;
								}
								TFloat distance = std::sqrt(std::max(static_cast<TFloat>(0), distance2));
								if (distance <= best_distance){
									best_distance = distance;
									this->get_L()[i] = ci;
								}
							}
							if (this->get_L()[i] != oldlabel){
								this->nchanges += 1;
							}
						}
						
						this->ndcalcs_X += this->ndata * this->ncentroids;
					}
				}
				
				
				
	

			};
		}
		
		virtual void set_C_tasks(){
			//set C
			this->C_tasks = {
				
				
		
				[this](TInt ti){
					if (ti == 0){
						std::vector<TFloat> simplesums (this->getncentroids()*this->getdimension(),0);
						std::vector<TInt> simplecounts (this->getncentroids(),0);
						for (TInt i = 0; i < this->getndata(); ++i){
							for (TInt di = 0; di < this->getdimension(); ++di){
								simplesums[this->get_L()[i]*this->getdimension() + di] += this->getdata()[i*this->getdimension() + di];
							}
							simplecounts[this->get_L()[i]]+=1;
						}
						for (TInt ci = 0; ci < this->getncentroids(); ++ ci){
							for (TInt di = 0; di < this->getdimension(); ++di){
								if (simplecounts[ci] != 0){
									this->get_C()[ci*this->getdimension() + di] = simplesums[ci*this->getdimension() + di] / static_cast<TFloat> (simplecounts[ci]);
								}
							}
						}
					}
				}
				
				


					
			};
		}
};

}

#endif
