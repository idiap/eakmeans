/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef PLL_HAMERLYKMEANS_13V0_H
#define PLL_HAMERLYKMEANS_13V0_H

#include "basehamerly.h"

//aka the Annulus algorithm
namespace kmeans{

template <typename TInt, typename TFloat>
void update_L_lower_upper_S_H_13v0(TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C,  const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const delta_C,  
/* The additions to base */
const TFloat * const data_l2s, const TFloat * const C_l2s, const TInt * const unordered_to_ordered, const TFloat * const C_ordered, const TFloat * const C_l2s_ordered, const TFloat * const C_l22s_ordered,
/* Additional to update */
TInt * const L2,
/* Usual to update */
 TInt * const L, TFloat * const lower, TFloat * const upper,  const TInt & round){
	
	nchanges = 0;
	ndcalcs = 0; 


	TFloat m;
	TInt oldlabel;
	std::unique_ptr<TFloat []> distances (new TFloat [ncentroids]);

	TFloat distance_to_L2;
		
	TInt insertion_lower;
	TInt insertion_upper;
	TInt nd_to_calc;
	TFloat searchradius;
	TInt L_ord;
	TInt L2_ord;
	
	//TODO: Hamerly checks that the label of data point is not max-mover (if fail, add second biggest budge). 
	TFloat max_deltaC_previous_round;
	TInt index_max_deltaC_previous_round;
	arrutilv2::set_argmaxmax(ncentroids, delta_C, index_max_deltaC_previous_round, max_deltaC_previous_round);
	
	for (TInt i = 0; i < ndata; ++i){
		lower[i] -= max_deltaC_previous_round;
		upper[i] += max_deltaC_previous_round;	
		m = std::max(halfminCC[L[i]], lower[i]);
		if (upper[i] > m){
			arrutilv2::set_l2(dimension, data + i*dimension,  C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upper[i], ndcalcs);		
			if (upper[i] > m){
				oldlabel = L[i];
				//the annulus trick goes here
				arrutilv2::set_l2(dimension, data + i*dimension,  C + L2[i]*dimension, data_l22s[i], C_l22s[L2[i]], distance_to_L2, ndcalcs);
				searchradius = std::max(upper[i], distance_to_L2) + 0.001;
				//std::cout << "Search Radius " << searchradius << std::endl;
				
				insertion_lower = std::max(static_cast<TInt>(0), arrutilv2::get_insertion_index_binarial(data_l2s[i] - searchradius, ncentroids, C_l2s_ordered));
				if (insertion_lower > 0){ insertion_lower -=1;}
				insertion_upper = std::max(static_cast<TInt>(0), arrutilv2::get_insertion_index_binarial(data_l2s[i] + searchradius, ncentroids, C_l2s_ordered) + 2);
	
				nd_to_calc = std::min(ncentroids, insertion_upper) - insertion_lower;
				//insertion_lower = 0;
				//nd_to_calc = ncentroids;
				arrutilv2::set_rl2s(dimension, data + i*dimension, nd_to_calc, C_ordered + insertion_lower*dimension, data_l22s[i], C_l22s_ordered + insertion_lower, distances.get(), ndcalcs);


				if (nd_to_calc >= 2){
					arrutilv2::set_argmin2min2nocheck(nd_to_calc, distances.get(), L_ord, L2_ord, upper[i], lower[i]);
				}
				else{
					
					std::cout << "error at index " << i << std::endl;
					std::cout << "nd to calc " << nd_to_calc  << std::endl;
					std::cout << "insertion " << insertion_lower << " " << nd_to_calc << " " << insertion_upper << std::endl;
					std::cout << "distance to nearest " << upper[i] << std::endl;
					std::cout << "distance to L2 " << distance_to_L2 << std::endl;
					std::cout << "search radius " << searchradius << std::endl;	
					
					std::cout << "data l2s[i] " << data_l2s[i] << std::endl;
					std::cout << "centroid l2s " << std::endl;
					for (TInt ci = 0; ci < ncentroids; ++ci){
						std::cout << C_l2s_ordered[ci] << " ";
					}
					std::cout << std::endl;
					
					
					throw std::logic_error("too few dists to calc, seems odd");
				}
				
				L[i] = unordered_to_ordered[insertion_lower + L_ord];
				L2[i] = unordered_to_ordered[insertion_lower + L2_ord];
				
				if (L[i] != oldlabel){
					++nchanges;
					++H[L[i]];
					--H[oldlabel];
					arrutilv2::addto(dimension, data + i*dimension, S + dimension*L[i]);
					arrutilv2::subtractfrom(dimension, data + i*dimension, S + dimension*oldlabel);
				}
			}
		}
	}
}

template <typename TInt, typename TFloat>
class P13V0 : public kmeans::BaseHamerly<TInt, TFloat>{
	
	private:
		std::unique_ptr<TFloat []> data_l2s;
		std::unique_ptr<TFloat []> C_l2s;
		std::unique_ptr<TInt []> unordered_to_ordered;
		std::unique_ptr<TFloat []> C_ordered;
		std::unique_ptr<TFloat []> C_l2s_ordered;
		std::unique_ptr<TFloat []> C_l22s_ordered;
		std::unique_ptr<TInt []> L2;
			
	protected:


		TFloat * const get_data_l2s(){
			return this->data_l2s.get();
		}
			
	
		TFloat * const get_C_l2s(){
			return this->C_l2s.get();
		}
	
		TInt * const get_unordered_to_ordered(){
			return this->unordered_to_ordered.get();
		}
		
		TFloat * const get_C_ordered(){
			return this->C_ordered.get();
		}
		
		TFloat * const get_C_l2s_ordered(){
			return this->C_l2s_ordered.get();
		}
		
		TFloat * const get_C_l22s_ordered(){
			return this->C_l22s_ordered.get();
		}
		
		TInt * const get_L2(){
			return L2.get();
		}
	
		TFloat * const get_lower(){
			return this->get_lower_base();
		}
		
		TFloat * const get_upper(){
			return this->get_upper_base();
		}
			

		
		std::function<void(TInt)> update_L_lower_upper_S_H_13v0_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				
				this->pll_principal_X(update_L_lower_upper_S_H_13v0<TInt, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_delta_C(), this->get_data_l2s() + x0, this->get_C_l2s(), this->get_unordered_to_ordered(), this->get_C_ordered(), this->get_C_l2s_ordered(), this->get_C_l22s_ordered(),  this->get_L2() + x0, this->get_L() + x0,  this->get_lower() + x0, this->get_upper() + x0, this->round);
			};
		}
		
	public:
		typedef kmeans::BaseHamerly<TInt, TFloat> BH;
		template<typename... Args>
		P13V0(Args&&... args): BH(std::forward<Args>(args)...),
		data_l2s { new TFloat [this->getndata()] },
		C_l2s { new TFloat [this->getncentroids()] },
		unordered_to_ordered { new TInt [this->getncentroids()] },
		C_ordered { new TFloat [this->getncentroids()*this->getdimension()] },
		C_l2s_ordered { new TFloat [this->getncentroids()] },
		C_l22s_ordered { new TFloat [this->getncentroids()] },
		L2 { new TInt [this->getndata()] }
		
		
		{
			this->setalgname("p13v0");
		}
		virtual ~P13V0(){}

		virtual void verbose_write_additional(){
			BH::verbose_write_additional();
			/* anything else to print ? */
		}
		
		std::function<void(TInt)> set_L2_ati(){ /* TODO : do what this function is supposed to : fill with second smallest. */
			return [this](TInt ti){
				if (ti == 0){
					for (TInt i = 0; i < this->getndata(); ++i){
						if (this->get_L()[i] == 0){
							this->get_L2()[i] = 1;
						}
						else{
							this->get_L2()[i] = 0;
						}
					}
				}	
			};
		}	

		virtual void set_initialisation_tasks(){
			
			
			this->initialisation_tasks = this->makeset_C_C_l22s_L_inds0_lower_upper_S_H_mati();
			
			/* set l2s, quick serial hack */
			this->initialisation_tasks.push_back(
				[this](TInt ti){
					if (ti == 0){
						for (TInt i = 0; i < this->getndata(); ++ i){
							this->data_l2s[i] = std::sqrt(std::max(static_cast<TFloat>(0), (this->get_data_l22s()[i])));
						}
					}
				}
			);
				
			this->initialisation_tasks.push_back(this->set_L2_ati());
		
					
			/* usual note about CC and halfminCC not needing to be set yet, no kmeans++ yet spiel bla de bla */
			
			}
		
		
		/*TODO : implement  */
		std::function<void(TInt)> update_unordered_to_ordered_C_C_l2_ati(){
			return [this](TInt ti){
				if (ti == 0){ //TODO: pllise this properly. 
					sort::set_argsorted_increasing(this->getncentroids(), this->get_C_l2s(), this->get_unordered_to_ordered());
					for (TInt ci = 0; ci < this->getncentroids(); ++ci){
						
						this->get_C_l22s_ordered()[ci] =  this->get_C_l22s()[this->unordered_to_ordered[ci]];
						this->get_C_l2s_ordered()[ci] = std::sqrt(std::max(static_cast<TFloat>(0), this->get_C_l22s_ordered()[ci]));
						
						std::memcpy(this->get_C_ordered() + ci*this->getdimension(), this->get_C() + this->unordered_to_ordered[ci]*this->getdimension(), sizeof(TFloat)*this->getdimension());
					}
					
					//arrutilv2::update_unordered_to_ordered_C_C_l22(this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_l22s(), this->get_unordered_to_ordered(), this->get_C_ordered(), this->get_C_l22s_ordered());
				}
			};
		}
			
		
		virtual void set_C_tasks(){
			this->C_tasks = {
				arrutilv2::update_C_C_l22s_delta_C_from_SH_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_sums(), this->get_counts(), this->get_C(), this->get_C_l22s(), this->get_delta_C(), this->ndcalcs_notX),

				arrutilv2::update_CC_halfminCC_ati(this->getnthreads(), this->getncentroids(), this->getdimension(), this->get_C(), this->get_C_l22s(), this->get_CC(), this->get_halfminCC(), this->ndcalcs_notX),
				
				/* set C l2s, quick serial hack */
				[this](TInt ti){
					if (ti == 0){
						for (TInt ci = 0; ci < this->getncentroids(); ++ ci){
							this->C_l2s[ci] = std::sqrt(std::max(static_cast<TFloat>(0), (this->get_C_l22s()[ci])));
						}
					}
				},
				
				
				
				this->update_unordered_to_ordered_C_C_l2_ati()
				
			};
		}
		
		virtual void set_X_tasks(){
					
			this->X_tasks = {
				/* TODO : implement the below function */
				this->update_L_lower_upper_S_H_13v0_ati()	
			};
		}
};

}

#endif



