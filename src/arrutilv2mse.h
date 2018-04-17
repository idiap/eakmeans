/*
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <james.newling@gmail.com>
All rights reserved.

eakmeans is a library for exact and approximate k-means written in C++ and
Python. This file is part of eakmeans. See file COPYING for more details.

This file is part of eakmeans.

eakmeans is free software: you can redistribute it and/or modify
it under the terms of the 3-Clause BSD Licence. See
https://opensource.org/licenses/BSD-3-Clause for more details.

eakmeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
COPYING for more details.
*/

#ifndef ARRUTILV2MSE_H
#define ARRUTILV2MSE_H

#include <cmath>
#include <memory>



namespace arrutilv2{

template <typename TInt, typename TSize, typename TFloat>
TFloat get_mse(TSize ndata, TInt ncentroids, const TFloat * const distances2, const TInt * const labels){
	TSize s_ncentroids = static_cast<TSize> (ncentroids);
	TFloat suml2s = 0;
	TFloat d2;
	for (TSize i = 0; i < ndata; ++i){
		d2 = distances2[i*s_ncentroids + labels[i]];
		if (d2 > 0){
			suml2s += std::sqrt(d2);
		}
		else if (d2 <-1e-5){
			throw std::runtime_error("negative value in get_mse of magnitude less than 1e-5. Probable cause: user has provided negative distance squared value");
		}
		else{
			//assumed rounding error (provide warning?)
		}
	}
	
	return suml2s/static_cast<TFloat> (ndata);
}

template <typename TSize, typename TFloat>
TFloat get_mse(TSize ndata, TFloat * distances){
		TFloat mse;
		set_row_sum_squares(static_cast<TSize>(1), ndata, distances, &mse);
		mse/=ndata;
		return mse;
}


//[ 0.5 * ridgeterm * sum (count - mean count)^2 ] / ndata
template <typename TSize, typename TFloat>
TFloat get_meanridge(TFloat ridgeterm, TSize ncentroids, TSize * const counts){
	
	auto ndata = 0;
	for (TSize si = 0; si < ncentroids; ++si){
		ndata += counts[si];
	}//std::accumulate(counts, counts + ncentroids);
	
	TFloat meancount = static_cast<TFloat> (ndata ) / static_cast<TFloat > (ncentroids);
	
	TFloat ridge_penalty = 0.;
	for (TSize ci = 0; ci < ncentroids; ++ ci){
		ridge_penalty += (counts[ci] - meancount)*(counts[ci] - meancount);
	}
	//TFloat ridge_penalty = static_cast<TFloat> (ridge_penalty_st);
	ridge_penalty *= ridgeterm/2.;
	ridge_penalty /= static_cast<TFloat> (ndata);
	return ridge_penalty;
}



template <typename TInt, typename TFloat>
TFloat getmeanl22at(TInt ncentroids, TInt dimension, const TFloat * const centroids, TInt ndata, const TFloat * const data, const TInt * const labels, const TFloat * const centroid_l22s, const TFloat * const data_l22s){		
	TFloat sum_variances = 0;
	for (TInt i = 0; i < ndata; ++ i){
		TFloat variance = 0;		
		for (TInt d = 0; d < dimension; ++ d){
			variance += data[i*dimension + d]*centroids[labels[i]*dimension + d];
		}
		variance *= -2;
		variance += data_l22s[i];
		variance += centroid_l22s[labels[i]];
		sum_variances += variance;
	}
	TFloat variance_estimate = sum_variances/static_cast<TFloat>(ndata);
	return variance_estimate;
}

//mse +  meanridgeerror
template <typename TSize, typename TFloat>
TFloat get_mse_ridge(TSize ndata, TFloat * distances, TFloat ridgeterm, TSize ncentroids, TSize * const counts){
	TFloat mse;
	set_row_sum_squares(static_cast<TSize>(1), ndata, distances, &mse);
	mse/=ndata;	
	TFloat meanridge = get_meanridge(ridgeterm, ncentroids, counts);	
	return mse + meanridge;
}





template <typename TInt, typename TFloat> 
TFloat get_sse_batchwise(TInt ndata, TInt nperbatch, TInt dimension, const TFloat * const data, TInt ncentroids, const TFloat * const centroids, const TFloat * const data_l22s, const TFloat * const centroid_l22s, TInt & ndcalcs){
	
	TInt nfullbatches = ndata/nperbatch;
	TInt nfinalbatch = ndata - nfullbatches*nperbatch;
	std::unique_ptr<TFloat []> distances_squared (new TFloat [nperbatch*ncentroids]);
	//data from the full batches
	TFloat sse = 0;
	for (TInt bi = 0; bi < nfullbatches; ++bi){
		set_rrl22ss(nperbatch, dimension, data + bi*dimension*nperbatch, ncentroids, centroids, data_l22s +bi*nperbatch, centroid_l22s, distances_squared.get());	
		for (TInt i = nperbatch*bi; i < nperbatch*(bi + 1); ++ i){
			sse += *std::min_element(distances_squared.get() + (i - nperbatch*bi)*ncentroids,  distances_squared.get() + (i - nperbatch*bi + 1)*ncentroids); 
		}
	}
	//data from the tail
	set_rrl22ss(nfinalbatch, dimension, data + nfullbatches*dimension*nperbatch, ncentroids, centroids, data_l22s + nfullbatches*nperbatch, centroid_l22s, distances_squared.get());	
	
	for (TInt i = nperbatch*nfullbatches; i < ndata; ++ i){
		sse += *std::min_element(distances_squared.get() + (i - nperbatch*nfullbatches)*ncentroids, distances_squared.get() + (i - nperbatch*nfullbatches + 1)*ncentroids); 
	}
	
	return sse;
}



		
}


#endif
