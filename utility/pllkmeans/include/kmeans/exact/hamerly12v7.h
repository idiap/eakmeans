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

#ifndef PLL_HAMERLYKMEANS_12V7_H
#define PLL_HAMERLYKMEANS_12V7_H

#include "baseexponion.h"
#include "sortutil.h"

namespace kmeans{



template <typename TInt, typename TTau, typename TFloat>
void update_L_lut_S_H_12v7(
TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C,  const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const delta_C, const TFloat * const max_deltaC_since,
TInt npartitions, std::pair<TFloat, TInt> * geometricpairs_halvies, TFloat * const partitionvalues_halvies, TInt * const geometricindices,
TInt * const L, TTau * const tau_lower, TFloat * const lower_at_last, TFloat * const upper,  const TTau & modround){

	if (npartitions != static_cast<TInt>(std::floor(std::log2(static_cast<double>(ncentroids - 1))))){
		throw std::runtime_error("npartitions is not as expected 12v7, throwing");
	}	
	nchanges = 0;
	ndcalcs = 0;
	TFloat m;	
	TInt base_index1;
	TFloat base_distance1;
	TFloat base_distance2;
	TInt insertion_index;
	TInt delta_part;
	TInt n_distances_to_calculate;
	TFloat lower;	
	TInt oldlabel;
	std::unique_ptr<TFloat []> distances (new TFloat [ncentroids]);
	for (TInt i = 0; i < ndata; ++i){
		lower = lower_at_last[i] - max_deltaC_since[tau_lower[i]];
		upper[i] += delta_C[L[i]];
		m = std::max(halfminCC[L[i]], lower);
		if (upper[i] > m){
			arrutilv2::set_l2(dimension, data + i*dimension,  C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upper[i], ndcalcs);		
			if (upper[i] > m){
				oldlabel = L[i];
				insertion_index = arrutilv2::get_insertion_index((upper[i] + halfminCC[L[i]]), npartitions - 1, partitionvalues_halvies + (npartitions - 1)*L[i]);
				tau_lower[i] = modround;
				if (insertion_index == npartitions - 1){//do as hamerly
					arrutilv2::set_rl2s(dimension, data + i*dimension, ncentroids, C, data_l22s[i], C_l22s, distances.get(), ndcalcs);
					arrutilv2::set_argminmin2nocheck(ncentroids, distances.get(), L[i], upper[i], lower_at_last[i]);
				}
				else{
					//don't remove this:
					if (insertion_index == 0){
						throw std::logic_error("(12v7) insertion index is 0, let's reflect on this for a moment. does this not mean that there are no centroids with radius upper + halfminCC[L[i]] of C[L[i]] ? which is a contradiction, as upper > halfminCC[L[i]] and so upper + halfminCC[L[i]] > minCC[L[i]], but by definition of minCC[L[i]] there is another  centroid within radius minCC[L[i]] of L[i]  ? This is a logic error, rethink code.");
					}
					//see 12v6 for what the following does. //TODO unify 12v6 and 12v7 codes
					delta_part = 1;
					n_distances_to_calculate = 0;
					for (TInt psi = 0; psi < insertion_index; ++psi){
						delta_part*=2;
						n_distances_to_calculate += delta_part;
					}
					distances[0] = upper[i];
					arrutilv2::set_l2s_at(n_distances_to_calculate, geometricindices + L[i]*(ncentroids - 1), dimension,  data + i*dimension, C, data_l22s[i], C_l22s, distances.get() + 1, ndcalcs);
					arrutilv2::set_argminmin2nocheck(n_distances_to_calculate + 1, distances.get(), base_index1, base_distance1, base_distance2);
					lower_at_last[i] = base_distance2;
					if (base_index1 != 0){
						L[i] = geometricindices[L[i]*(ncentroids - 1) + base_index1 -1];
						upper[i] = base_distance1;
					}
					/* no change of label */
					else{
					}					
				}
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








	
template <typename TInt, typename TTau, typename TFloat>
class P12V7 : public kmeans::BaseExponion<TInt, TFloat>{
	
	private:

		TTau cutperiod;		
		std::unique_ptr<TTau []> tau_lower;
		std::unique_ptr<TFloat []>  cumabs;
		std::unique_ptr<TFloat []>  max_deltaC_since;
		
		//(1) done always (cutrounds and non-cutrounds)
		std::function<void(TInt)> update_cumabs_max_ati(){
			return [this](TInt ti){
				TInt modround = 1 + (this->round  - 1) % this->get_cutperiod();
				TInt r0 = (ti*modround)/this->getnthreads();
				TInt r1 = ((ti + 1)*modround)/this->getnthreads();
				arrutilv2::update_cumabs_max(this->getncentroids(), this->get_delta_C(), r1 - r0, this->get_cumabs() + r0*this->getncentroids(), this->get_max_deltaC_since() + r0);
					
			};		
		}
		

		void reset_tau_lowers_ati(TInt ti){
			TInt x0 = (ti*this->getndata())/this->getnthreads();
			TInt x1 = ((ti+1)*this->getndata())/this->getnthreads();
			for (TInt i = x0; i < x1; ++i){
				this->get_lower_at_last()[i] -= this->max_deltaC_since[this->tau_lower[i]];
				this->tau_lower[i] = 0;
			}
		}

		//(2) if in cutround 		
		std::function<void(TInt)> reset_tau_lowers_if_cut_ati(){
			return [this](TInt ti){
				if (this->round%this->cutperiod == 0){
					reset_tau_lowers_ati(ti);
				}
			};
		}
		
		
		void reset_cumabs_max_ati(TInt ti){
			TInt r0 = (ti*(this->get_cutperiod() + 1))/this->getnthreads();
			TInt r1 = ((ti + 1)*(this->get_cutperiod() + 1))/this->getnthreads();
			for (TInt r = r0; r < r1; ++r){
				std::fill_n(this->get_cumabs() + r0*this->getncentroids(), (r1 -r0)*this->getncentroids(),0 );
				std::fill_n(this->get_max_deltaC_since() + r0, (r1 -r0), 0);
			}
		}
		
		//(3) if in cutround 		
		std::function<void(TInt)> reset_cumabs_max_if_cut_ati(){
			return [this](TInt ti){
				if (this->round%this->cutperiod == 0){
					reset_cumabs_max_ati(ti);
				}
			};
		}
		
			
	protected:
		
		TTau get_cutperiod(){
			return cutperiod;
		} 
		
		TTau * const get_tau_lower(){
			return tau_lower.get();
		}
		
		TFloat * const get_cumabs(){
			return cumabs.get();
		}
		
		TFloat * const get_max_deltaC_since(){
			return max_deltaC_since.get();
		}
		
		TFloat * const get_lower_at_last(){
			return this->get_lower_base();
		}
		
		TFloat * const get_upper(){
			return this->get_upper_base();
		}
				
			
		std::function<void(TInt)> update_L_lut_S_H_12v7_ati(){
			return [this](TInt ti){
				
				TInt modround = this->round % this->get_cutperiod();		
				
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				
				this->pll_principal_X(kmeans::update_L_lut_S_H_12v7<TInt, TTau, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_delta_C(), this->get_max_deltaC_since(), this->get_npartitions(), this->get_geometricpairs_halvies(), this->get_partitionvalues_halvies(), this->get_geometricindices(), this->get_L() + x0,  this->get_tau_lower() + x0, this->get_lower_at_last() + x0, this->get_upper() + x0, modround);
				
			};
		}

		
		
	public:
		typedef kmeans::BaseExponion<TInt, TFloat> BEX;
		template<typename... Args>
		P12V7(Args&&... args): BEX(std::forward<Args>(args)...),
		
	/* ignoring constants, M = C^2 + Nd + cC
	 * where M is memory use, N is ndata, d is dimension, C in ncentroids, c is cutperiod
	 * to ensure that memory is no larger (O) than unsmartbound exponion (12v6),  
	 * cC < max (C^2, Nd) --> c < max (C, Nd/C)
	 * we add the further bound that 1 <= c <= 50, resulting in the expression which follows.
	 * */
		cutperiod { static_cast<TTau> (
		std::max(static_cast<TInt>(1), ( //step 3 : make greater than or equal to 1
		std::min(static_cast<TInt>(50),  //step 2: make less than or equal to 50
		std::max(this->getncentroids(), static_cast<TInt>(
		static_cast<TFloat>(this->getndata()*this->getdimension())/static_cast<TFloat>(this->getncentroids()) //step 1: set to max (C, Nd/C)
		)))))) },

		tau_lower { new TTau [this->getndata()] }, 
		//not sure whether cuptperiod, cutperiod + 1 etc.
		cumabs {new TFloat [this->getncentroids()*(cutperiod + 1)]},
		max_deltaC_since {new TFloat [this->getncentroids()*(cutperiod + 1)]}

		{
			std::fill_n(tau_lower.get(), this->getndata(), 0);
			std::fill_n(cumabs.get(), this->getncentroids()*(cutperiod + 1), 0);
			std::fill_n(max_deltaC_since.get(), (cutperiod + 1), 0);
			//if (this->get_cout_verbosity() > 0){
				//std::cout << "cutperiod set to " << cutperiod << std::endl;
			//}
			this->setalgname("p12v7");
		}
			
		virtual ~P12V7(){}

		virtual void verbose_write_additional(){
		BEX::verbose_write_additional();
			/* anything else to print ? */
		}

		virtual void set_initialisation_tasks(){
			this->exponion_set_initialisation_tasks();
		}
		
		virtual void set_C_tasks(){
			this->exponion_set_C_tasks();
			//for updating cumabs and max_delta_C_since
			this->C_tasks.push_back (this->update_cumabs_max_ati());
			//the following is not really a C task, but has to be here for chronological consistency (C_tasks and X_tasks just get concatenated, syntactic gloss really)
			this->C_tasks.push_back (this->reset_tau_lowers_if_cut_ati());  
							
			this->C_tasks.push_back (this->reset_cumabs_max_if_cut_ati());
			
			/* barriers 'do not' slow things down , the following does not affect runtime : 
			
			for (TInt x = 0; x < 13; ++x){
				this->C_tasks.push_back ([](TInt ti){});
			}
			* 
			* */
		}
		
		
		virtual void set_X_tasks(){
			this->X_tasks = {
				this->update_L_lut_S_H_12v7_ati()
			};
		}
};

}

#endif

