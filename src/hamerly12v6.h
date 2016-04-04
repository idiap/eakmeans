#ifndef PLL_HAMERLYKMEANS_12V6_H
#define PLL_HAMERLYKMEANS_12V6_H

#include "baseexponion.h"
#include "sortutil.h"

namespace kmeans{



template <typename TInt, typename TFloat>
void update_L_lower_upper_S_H_12v6(
TInt ncentroids, TInt dimension, TFloat * const S, TInt * const H , TInt & nchanges, TInt &ndcalcs, 
TInt ndata, const TFloat * const data, const TFloat * const C,  const TFloat * const data_l22s, const TFloat * const C_l22s, const TFloat * const CC, const TFloat * const halfminCC, const TFloat * const delta_C, 
TInt npartitions, std::pair<TFloat, TInt> * geometricpairs_halvies, TFloat * const partitionvalues_halvies, TInt * const geometricindices,
TInt * const L, TFloat * const lower, TFloat * const upper,  const TInt & round){
		
	if (npartitions != static_cast<TInt>(std::floor(std::log2(static_cast<double>(ncentroids - 1))))){
		throw std::runtime_error("npartitions is not as expected, what's up? If being called by 12v6 class, logic error. If user called, it may be a runtime_error");
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
	TInt oldlabel;	
		
	
	std::unique_ptr<TFloat []> distances (new TFloat [ncentroids]);
	TFloat max_deltaC_previous_round = *std::max_element(delta_C, delta_C + ncentroids);
	
	for (TInt i = 0; i < ndata; ++i){
		lower[i] -= max_deltaC_previous_round; //lower bound on second nearest
		upper[i] += delta_C[L[i]]; //upper bound on nearest (L[i])
		m = std::max(halfminCC[L[i]], lower[i]);
		if (upper[i] > m){
			arrutilv2::set_l2(dimension, data + i*dimension,  C + L[i]*dimension, data_l22s[i], C_l22s[L[i]], upper[i], ndcalcs);		
			if (upper[i] > m){
				oldlabel = L[i];
				//why is this correct ? see the paper, and note partionvalues_HALVIES :)
				insertion_index = arrutilv2::get_insertion_index((upper[i] + halfminCC[L[i]]), npartitions - 1, partitionvalues_halvies + (npartitions - 1)*L[i]); // making seach binarial does not help 
				//insertion_index = arrutilv2::get_insertion_index((upper[i] + m), npartitions - 1, partitionvalues_halvies + (npartitions - 1)*L[i]);
				if (insertion_index == npartitions - 1){//do as hamerly
					arrutilv2::set_rl2s(dimension, data + i*dimension, ncentroids, C, data_l22s[i], C_l22s, distances.get(), ndcalcs);
					arrutilv2::set_argminmin2nocheck(ncentroids, distances.get(), L[i], upper[i], lower[i]);
				}
				else{
					//TODO: If i remove this check, runs faster, bizarre
					if (insertion_index == 0){
						throw std::logic_error("insertion index is 0, let's reflect on this for a moment. does this not mean that there are no centroids with radius upper + halfminCC[L[i]] of C[L[i]] ? which is a contradiction, as upper > halfminCC[L[i]] and so upper + halfminCC[L[i]] > minCC[L[i]], but by definition of minCC[L[i]] there is another  centroid within radius minCC[L[i]] of L[i]  ? This is a logic error, rethink code.");
					}
					/* number of centroids to whom distance needs to be calculated for different values of insertion_index
					 * 0 : error
					 * 1 : 2  : (nearest and 2 nd nearest)
					 * 2 : 6  : (nearest -> 6  th nearest)
					 * 3 : 14 : (nearest -> 14 th nearest)
					 * this is calculated in the following for loop
					 */
					delta_part = 1;
					n_distances_to_calculate = 0;
					for (TInt psi = 0; psi < insertion_index; ++psi){
						delta_part*=2;
						n_distances_to_calculate += delta_part;
					}
					distances[0] = upper[i];
					arrutilv2::set_l2s_at(
					n_distances_to_calculate, 
					geometricindices + L[i]*(ncentroids - 1), 
					dimension,  
					data + i*dimension, 
					C,  
					data_l22s[i], 
					C_l22s, distances.get() + 1, ndcalcs);
					
					arrutilv2::set_argminmin2nocheck(n_distances_to_calculate + 1, distances.get(), base_index1, base_distance1, base_distance2);					
					if (base_index1 != 0){
						L[i] = geometricindices[L[i]*(ncentroids - 1) + base_index1 -1];
						upper[i] = base_distance1;
						lower[i] = base_distance2;
					}
					/* no change of label */
					else{
						lower[i] = base_distance2;
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
	
template <typename TInt, typename TFloat>
class P12V6 : public kmeans::BaseExponion<TInt, TFloat>{
	
	private:
			
	protected:
		
		TFloat * const get_lower(){
			return this->get_lower_base();
		}
		
		TFloat * const get_upper(){
			return this->get_upper_base();
		}
			
		std::function<void(TInt)> update_L_lower_upper_S_H_12v6_ati(){
			return [this](TInt ti){
				TInt x0 = (ti*this->getndata())/this->getnthreads();
				this->pll_principal_X(kmeans::update_L_lower_upper_S_H_12v6<TInt, TFloat>, ti, this->get_CC(), this->get_halfminCC(), this->get_delta_C(), this->get_npartitions(), this->get_geometricpairs_halvies(), this->get_partitionvalues_halvies(), this->get_geometricindices(), this->get_L() + x0,  this->get_lower() + x0, this->get_upper() + x0, this->round);
			};
		}

		
		
	public:
		typedef kmeans::BaseExponion<TInt, TFloat> BEX;
		template<typename... Args>
		P12V6(Args&&... args): BEX(std::forward<Args>(args)...)

		{
			this->setalgname("p12v6");
		}
		virtual ~P12V6(){}

		virtual void verbose_write_additional(){
			BEX::verbose_write_additional();
			/* anything else to print ? */
		}

		virtual void set_initialisation_tasks(){
			this->exponion_set_initialisation_tasks();
		}
		
		virtual void set_C_tasks(){
			this->exponion_set_C_tasks();
		}
		
		virtual void set_X_tasks(){
			this->X_tasks = {
			this->update_L_lower_upper_S_H_12v6_ati()	
			};
		}
};

}

#endif

