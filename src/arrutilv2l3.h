#ifndef ARRUTILV2L3_H
#define ARRUTILV2L3_H

/* for c++11 threading */

#include "arrutilv2l0.h"
#include "arrutilv2l1.h"
#include "arrutilv2l2.h"
#include "arrutilv2copy.h"
#include "arrutilv2discrete.h"

#include <mutex>
#include <thread>

namespace arrutilv2{
	
	
/* ati functions return functions of form std::function<void(TInt)>, returned item being bound to data, and acting at TInt index
 * 
 * for updates where there is no mutex needed for parallel versions, the function layering goes like:
 *  - fn_update (in l2)
 *  - - fn_update_ati (in l3, returns version of l2 fn_update bound to data pointers) 
 * 
 * for updates where there is a mutex needed for parallel versions, the function layering goes like this:
 *  - fn_update (in l2)
 *  - - fn_update (in l3, calls l2 version and uses mutex)
 *  - - - fn_update_ati (in l3, returns version of l3 fn_update bound to data pointers)
 * 
 * 
 * see also mati functions in clusteringutil/initv2
 * 
 * */

template <typename TInt, typename TFloat>
/* update sums and counts according to change calculated above */	
/* this serialisation is not expensive : with 400 centroids 20000 data, repeating this locked step 10 times increasing running time ~2%. */	
inline void update_sums_counts(const TInt & ncentroids, const TInt & dimension, TFloat * const sums, const TFloat * const dsums, TInt * const counts, const TInt * const dcounts){
	
	for (TInt j = 0; j < ncentroids*dimension; ++j){
		sums[j] += dsums[j];
	}
	for (TInt j = 0; j < ncentroids; ++j){
		counts[j] += dcounts[j];
	}
}

	
template <typename TInt, typename TFloat>  
std::function<void(TInt)> update_C_C_l22s_from_SH_ati(TInt nthreads, TInt ncentroids, TInt dimension, const TFloat * const sums, const TInt * const counts, TFloat * const C, TFloat * const C_l22s){

	
	return [dimension, ncentroids, nthreads, C, C_l22s, sums, counts](TInt ti){
		TInt c0 = (ti*ncentroids)/nthreads;
		TInt c1 = ((ti+1)*ncentroids)/nthreads; 
		arrutilv2::update_C_C_l22s_from_SH(c1 - c0, dimension, sums + c0*dimension, counts + c0, C + c0*dimension, C_l22s + c0);};
}

template <typename TInt, typename TNum>
std::function<void(TInt)> set_vector_sum_ati(TInt nthreads, TInt dimension, const TNum * const A, const TNum * const B, TNum * const C){
	
		return [nthreads, dimension, A, B, C](TInt ti){
			TInt i0 = (ti*dimension)/nthreads;
			TInt i1 = ((ti+1)*dimension)/nthreads;
			std::memcpy(C + i0, A + i0, sizeof(TNum)*(i1 - i0));
			arrutilv2::addto(i1 - i0, B + i0, C + i0);  
		};
}


template <typename TInt, typename TFloat>  
std::function<void(TInt)> update_C_C_l22s_delta_C_from_SH_ati(TInt nthreads, TInt ncentroids, TInt dimension, const TFloat * const sums, const TInt * const counts, TFloat * const C, TFloat * const C_l22s, TFloat * const delta_C, std::atomic<TInt> & ndcalcs){

	return [nthreads, ncentroids, dimension, sums, counts, C, C_l22s, delta_C, &ndcalcs](TInt ti){
		TInt c0 = (ti*ncentroids)/nthreads;
		TInt c1 = ((ti+1)*ncentroids)/nthreads; 
		TInt local_ndcalcs = 0;
		update_C_C_l22s_delta_C_from_SH(c1 - c0, dimension, sums + c0*dimension, counts + c0, C + c0*dimension, C_l22s + c0, delta_C + c0, local_ndcalcs);
				
		ndcalcs += local_ndcalcs;
	};
}



template <typename TInt, typename TFloat>  
std::function<void(TInt)> update_delta_G_from_delta_C_ati(TInt ncentroids, const TFloat * const delta_C, TInt ngroups, const TInt * const groupparts, const TInt * const groupsizes, TFloat * const delta_G){
	
	/*TODO : split task, by groups */
	return [ncentroids, delta_C, ngroups, groupparts, delta_G](TInt ti){
		if (ti == 0){
			update_delta_G_from_delta_C(ncentroids, delta_C, ngroups, groupparts, delta_G);
		}
	};
	
}



template <typename TInt, typename TFloat>  
std::function<void(TInt)> update_C_C_hist_from_SH_ati(TInt nthreads, TInt ncentroids, TInt dimension, const TFloat * const sums, const TInt * const counts, TFloat * const C, TFloat * const C_l22s, TFloat * const C_hist, TFloat * const C_l22s_hist, TInt & round, TInt cutperiod){
	
	
	/* plonk C in hist and update C
	 * if cutperiod = 4, plonk positions
	 * (by round vertical) look like this,
	 * 
	 *   0 1 2 3
	 * 1 * - - -          
	 * 2 - * - - 
	 * 3 - - * - 
	 * 4 - - - *
	 * 5 * - - - 
	 * 6 - * - - 
	 * 7 - - * - 
	 * etc
	 * */
	return [nthreads, ncentroids, dimension, sums, counts, C, C_l22s, C_hist, C_l22s_hist, &round, cutperiod](TInt ti){
		TInt c0 = (ti*ncentroids)/nthreads;
		TInt c1 = ((ti+1)*ncentroids)/nthreads;
		TInt C_hist_s = ((round - 1)%cutperiod)*ncentroids*dimension; 
		TInt C_l22s_hist_s = ((round - 1)%cutperiod)*ncentroids;
		move_then_update_C_C_l22s_from_SH(c1 - c0, dimension, sums + c0*dimension, counts + c0, C + c0*dimension, C_l22s + c0, C_hist + C_hist_s + c0*dimension, C_l22s_hist + C_l22s_hist_s + c0);
	};
}





template <typename TInt, typename TFloat>  
std::function<void(TInt)> update_u_deltaC_from_C_C_hist_ati(TInt nthreads, TInt ncentroids, TInt dimension, const TFloat * const C,  const TFloat * const C_hist, TFloat * const u_deltaC,  TInt & round, TInt cutperiod,  std::atomic<TInt> & ndcalcs){

	/* update u_deltaC from C and hist
	 * if cutperiod = 4, the updates affectuated are
	 * 
	 *   0 1 2 3 4
	 * 1 * 0 - - -         
	 * 2 * * 0 - -
	 * 3 * * * 0 -
	 * 4 * * * * 0
	 * 5 * 0 - - -
	 * 6 * * 0 - -
	 * 7 * * * 0 -
	 * etc
	 * */
	 
		/* take C and C_hist and set u_deltaC to be the displacement between. 
		 * set element |C_hist| if u_deltaC to be zeros (corresponding to no move size history = present) 
		 * 
		 * when round%cutperiod == 1, history contains 1 centroid
		 * if round%cutperiod == 0, history has cutperiod centroids.
		 * */
				
	return [nthreads, ncentroids, dimension, C, C_hist, u_deltaC, &round, cutperiod, &ndcalcs](TInt ti){

		TInt nrounds_to_proc = 1 + (round - 1)%cutperiod;
		TInt r0 = (ti*nrounds_to_proc)/nthreads;
		TInt r1 = ((ti + 1)*nrounds_to_proc)/nthreads;
		
		TInt local_ndcalcs = 0;
		if (ti == nthreads - 1){
			update_u_deltaC_from_C_C_hist(ncentroids, dimension, C, r1 - r0, C_hist + r0*ncentroids*dimension, u_deltaC + r0*ncentroids, local_ndcalcs);
			std::fill_n(u_deltaC + r1*ncentroids, ncentroids, static_cast<TInt>(0));
		}
		else{
			update_u_deltaC_from_C_C_hist(ncentroids, dimension, C, r1 - r0, C_hist + r0*ncentroids*dimension, u_deltaC + r0*ncentroids, local_ndcalcs);
		}
		
		ndcalcs += local_ndcalcs;

	};
}



template <typename TInt, typename TFloat> 
void set_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, std::mutex & mut){

	
	arrutilv2::set_S_H(ndata, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dsums, dcounts);
	std::lock_guard<std::mutex> gluk(mut);
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);
}


template <typename TInt, typename TFloat>
std::function<void(TInt)> set_S_H_ati(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, const TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, std::mutex & work_mutex){

	return [dimension, ndata, ncentroids, nthreads,  data, C, data_l22s, C_l22s, L, dsums, dcounts, sums, counts, &work_mutex](TInt ti){
		TInt x0 = (ti*ndata)/nthreads;
		TInt x1 = ((ti+1)*ndata)/nthreads;	
		arrutilv2::set_S_H(x1-x0, dimension,  data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, L + x0, dsums + ti*dimension*ncentroids, dcounts + ti*ncentroids,  sums, counts, work_mutex);};
}


template <typename TInt, typename TFloat> 
void set_L_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, std::mutex & mut){
	arrutilv2::set_L_S_H(ndata, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dsums, dcounts);
	std::lock_guard<std::mutex> gluk(mut);
	update_sums_counts(mut, ncentroids, dimension, sums, dsums, counts, dcounts);
}






template <typename TInt, typename TFloat>
void update_L_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & mut){
	/* update labels and get change to sums and counts */
	TInt local_nchanges = 0;
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 


	
	arrutilv2::update_L_S_H(ndata, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dsums, dcounts, local_nchanges);

	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);
}
	



template <typename TInt, typename TFloat>
std::function<void(TInt)> update_L_S_H_ati(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & work_mutex, std::atomic<TInt> & ndcalcs){
	
	
	return [dimension, ndata, ncentroids, nthreads,  data, C, data_l22s, C_l22s, L, dsums, dcounts, sums, counts, &nchanges, &work_mutex, &ndcalcs](TInt ti){
		
		TInt x0 = (ti*ndata)/nthreads;
		TInt x1 = ((ti+1)*ndata)/nthreads;	
		arrutilv2::update_L_S_H(x1-x0, dimension,  data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, L + x0, dsums + ti*dimension*ncentroids, dcounts + ti*ncentroids,  sums, counts, nchanges, work_mutex);
		ndcalcs += (x1 - x0)*ncentroids;
		};
}


template <typename TInt, typename TFloat>
inline void update_L_dn_S_H_batch(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dn, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & mut){
	/* update labels and get change to sums and counts */
	TInt local_nchanges = 0;
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 
	
	arrutilv2::update_L_dn_S_H_batch(ndata, nperbatch, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dn, dsums, dcounts, local_nchanges);

	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);
}






template <typename TInt, typename TFloat>
void update_L_S_H_batch(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & mut){
	
	//The following appears to be as fast as using the function update_L_S_H_batch of level 2 (as opposed to the version which does dn as well) FALSE STATEMENT. 
	//std::unique_ptr<TFloat []>  dn (new TFloat [ndata]);
	//update_L_dn_S_H_batch(ndata, nperbatch, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dn.get(), dsums, dcounts, sums, counts, nchanges, mut);
	
//}


	/* update labels and get change to sums and counts */
	TInt local_nchanges = 0;
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 

	//TODO: replace with (uncomment, comment) (faster than following)
	arrutilv2::update_L_S_H_batch(ndata, nperbatch, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dsums, dcounts, local_nchanges);

	
	//std::unique_ptr<TFloat []>  dn (new TFloat [ndata]);
	//arrutilv2::update_L_dn_S_H_batch(ndata, nperbatch, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dn.get(), dsums, dcounts, local_nchanges);
	

	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);

}



template <typename TInt, typename TFloat>
void update_L_S_H_batch_increment_only(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & mut){
	/* update labels and get change to sums and counts */
	TInt local_nchanges = 0;
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 
	
	arrutilv2::update_L_S_H_batch_increment_only(ndata, nperbatch, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dsums, dcounts, local_nchanges);

	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);
}



//base on update_L_dn_S_H_batch_increment_only
template <typename TInt, typename TFloat> 
inline void set_L_lowers_dn_and_increment_S_H(TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const lowers, TFloat * const dn, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & mut){

	TInt local_nchanges = 0;
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 
	           
	arrutilv2::set_L_lowers_dn_and_increment_S_H(ndata, dimension, data, ncentroids, C, data_l22s, C_l22s, L, lowers, dn, dsums, dcounts, local_nchanges);

	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);	
}

template <typename TInt, typename TFloat>
void update_L_dn_S_H_batch_increment_only(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dn, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & mut){
	/* update labels and get change to sums and counts */
	TInt local_nchanges = 0;
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 
	
	arrutilv2::update_L_dn_S_H_batch_increment_only(ndata, nperbatch, dimension, data, ncentroids, C, data_l22s, C_l22s, L, dn, dsums, dcounts, local_nchanges);

	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);
}





template <typename TInt, typename TFloat>
std::function<void(TInt)> update_L_S_H_batch_ati(TInt nthreads, TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & work_mutex, std::atomic<TInt> & ndcalcs){	
	return [dimension, ndata, nperbatch, ncentroids, nthreads,  data, C, data_l22s, C_l22s, L, dsums, dcounts, sums, counts, &nchanges, &work_mutex, &ndcalcs](TInt ti){
		TInt x0 = (ti*ndata)/nthreads;
		TInt x1 = ((ti+1)*ndata)/nthreads;	
		
		arrutilv2::update_L_S_H_batch(x1-x0, nperbatch, dimension, data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, L + x0, dsums + ti*dimension*ncentroids, dcounts + ti*ncentroids,  sums, counts, nchanges, work_mutex);
		
		ndcalcs += (x1 - x0)*ncentroids;
		};
}


/* probs deja vu  */
template <typename TInt, typename TFloat>
std::function<void(TInt)> set_L_S_H_batch_ati(TInt nthreads, TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & work_mutex, std::atomic<TInt> & ndcalcs){	
	return [dimension, ndata, nperbatch, ncentroids, nthreads,  data, C, data_l22s, C_l22s, L, dsums, dcounts, sums, counts, &nchanges, &work_mutex, &ndcalcs](TInt ti){
		TInt x0 = (ti*ndata)/nthreads;
		TInt x1 = ((ti+1)*ndata)/nthreads;	
		arrutilv2::update_L_S_H_batch_increment_only(x1-x0, nperbatch, dimension, data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, L + x0, dsums + ti*dimension*ncentroids, dcounts + ti*ncentroids,  sums, counts, nchanges, work_mutex);
		ndcalcs += (x1 - x0)*ncentroids;
		};
}


////a switch on two preceding functions, used in GrowBatchMse. if an_update is true, performs update_L_S_H_batch, else performs update_L_S_H_batch_increment_only. 
//template <typename TInt, typename TFloat>
//std::function<void(TInt)> update_or_setandincr_L_S_H_batch_ati(TInt nthreads, TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt * const L, TFloat * const dsums, TInt * const dcounts, TFloat * const sums, TInt * const counts, TInt & nchanges, std::mutex & work_mutex, std::atomic<TInt> & ndcalcs, bool & an_update){
	//return [dimension, ndata, nperbatch, ncentroids, nthreads,  data, C, data_l22s, C_l22s, L, dsums, dcounts, sums, counts, &nchanges, &work_mutex, &ndcalcs](TInt ti){
		//TInt x0 = (ti*ndata)/nthreads;
		//TInt x1 = ((ti+1)*ndata)/nthreads;	
		//if (an_update == true){
			//arrutilv2::update_L_S_H_batch(x1-x0, nperbatch, dimension, data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, L + x0, dsums + ti*dimension*ncentroids, dcounts + ti*ncentroids,  sums, counts, nchanges, work_mutex);
		//}
		//else{
			//arrutilv2::update_L_S_H_batch_increment_only(x1-x0, nperbatch, dimension, data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, L + x0, dsums + ti*dimension*ncentroids, dcounts + ti*ncentroids,  sums, counts, nchanges, work_mutex);
		//}	
		//ndcalcs += (x1 - x0)*ncentroids;
	//};
//}





template <typename TInt, typename TFloat, typename Function, typename... Args>
void pll_update_L_etc(
/* a fuction whose first parameters are... */
const Function & update_L_etc, 
/* these guys, and whose subsequent parameters are... */
TInt ncentroids, TInt dimension, TFloat * const sums, TFloat * const dsums, TInt * const counts, TInt * const dcounts, TInt & nchanges, std::atomic<TInt> & ndcalcs, std::mutex & mut, 
/* these guys. */
Args&&... args
/* the function is a label (and bounds etc) updater, where whenever a label changes an adjustment needs to be made to global sums and counts. but to prevent races the global sums and counts need to be updated via a mutex.  */
)

{
	/* local n changes and nd calcs */
	TInt local_nchanges = 0;
	TInt local_ndcalcs = 0;
	
	/* fill local sums and counts with 0s */
	std::fill_n(dsums, ncentroids*dimension, 0);
	std::fill_n(dcounts, ncentroids, 0); 
	
	/* main function call */
	update_L_etc(ncentroids, dimension, dsums, dcounts, local_nchanges, local_ndcalcs, std::forward<Args>(args)...);	

	/* update globals (not that std::atomic<TInt> ndcalcs does not need to be under lock below. Before I had ndcalcs as a TInt, and was having race conditions arising, it was in an attempt to stop these races that I have it underlock below. In the end making ndcalcs std::atomic<TInt> was the cleanest solution though ) */
	std::lock_guard<std::mutex> gluk(mut);
	nchanges += local_nchanges;
	ndcalcs += local_ndcalcs;
	update_sums_counts(ncentroids, dimension, sums, dsums, counts, dcounts);
} 




template <typename TInt, typename TFloat>
std::function<void(TInt)> update_CC_halfminCC_ati(TInt nthreads, TInt ncentroids, TInt dimension, const TFloat * const C, const TFloat * const C_l22s,  TFloat * const CC, TFloat * const halfminCC, std::atomic<TInt> & ndcalcs){
	
	return [nthreads, ncentroids, dimension, C, C_l22s, CC, halfminCC, &ndcalcs](TInt ti){
		TInt c0 = (ti*ncentroids)/nthreads;
		TInt c1 = ((ti+1)*ncentroids)/nthreads;
		/* correct indices of CC and halfminCC calculated in partial function, not nec to increment pointer here */
		TInt local_ndcalcs = 0;
		update_CC_halfminCC_partial(c0, c1, ncentroids, dimension, C, C_l22s, CC, halfminCC, local_ndcalcs);
		ndcalcs += local_ndcalcs;
	};
}

template <typename TInt, typename TFloat>
std::function<void(TInt)> update_CC_ati(TInt nthreads, TInt ncentroids, TInt dimension, const TFloat * const C, const TFloat * const C_l22s,  TFloat * const CC, std::atomic<TInt> & ndcalcs){
	
	return [nthreads, ncentroids, dimension, C, C_l22s, CC, &ndcalcs](TInt ti){
		TInt c0 = (ti*ncentroids)/nthreads;
		TInt c1 = ((ti+1)*ncentroids)/nthreads;
		/* correct indices of CC and halfminCC calculated in partial function, not nec to increment pointer here */
		TInt local_ndcalcs = 0;		
		set_l2gramm_partial(c0, c1, ncentroids, dimension, C, C_l22s, CC, local_ndcalcs);
		ndcalcs += local_ndcalcs;
	};
}





template <typename TInt, typename TFloat>
std::function<void(TInt)> set_L_group_dn_glowers_ati(TInt nthreads, TInt ndata, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const C, const TFloat * const data_l22s, const TFloat * const C_l22s, TInt ngroups, const TInt * const groupparts, const TInt * const  groupsizes, TInt * const L, TInt * const group, TFloat * const dn,  TFloat * const glowers, std::atomic<TInt> & ndcalcs){

		
	return [nthreads, ndata, dimension, data, ncentroids, C, data_l22s, C_l22s, ngroups, groupparts, groupsizes, group, L , glowers, dn, &ndcalcs](TInt ti)		
	{
		TInt local_ndcalcs = 0;
		TInt x0 = (ti*ndata)/nthreads;
		TInt x1 = ((ti+1)*ndata)/nthreads; 
		arrutilv2::set_L_group_dn_glowers(x1 - x0, dimension, data + x0*dimension, ncentroids, C, data_l22s + x0, C_l22s, ngroups, groupparts, groupsizes, L + x0, group + x0, dn + x0, glowers + x0*ngroups, local_ndcalcs);
		ndcalcs += local_ndcalcs;
	};


}



		
		

template <typename TInt, typename TFloat>
std::function<void(TInt)> set_rl22s_ati(TInt nthreads, TInt ndata, TInt dimension,  const TFloat * const A, const TFloat * const B, TFloat * const l22s, std::atomic<TInt> & ndcalcs){
	
	return [nthreads, ndata, dimension, A, B, l22s, &ndcalcs](TInt ti){
		
		TInt local_ndcalcs = 0;
		TInt x0 = (ti*ndata)/nthreads;
		TInt x1 = ((ti+1)*ndata)/nthreads; 
		
		
		set_rl22s(x1 - x0, dimension, A + x0*dimension, B + x0*dimension, l22s + x0, local_ndcalcs);
		ndcalcs += local_ndcalcs;
		
	};
}



}

#endif






