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

#ifndef MINIBATCHAPP_H
#define MINIBATCHAPP_H

namespace minibatchapp{

template <typename TInt>
class MiniBatchApp{
	public:
		//nuber of data used for each centroid update
		TInt batchsize;
		//number of centroid updates in one complete round of data
		TInt nsubrounds;
		//this round mod nsubrounds
		TInt subround;
		//the amount of data used in the nsubround'th update, the other updates use batchsize datapoints
		TInt lastbatchsize;
		//The size of the first batch. 
		TInt initialising_batch_size;
		//Number of changes on each batch since previous time that batch processed. 
		std::vector<TInt> nchanges_on_batch;		


		MiniBatchApp() = default;

		MiniBatchApp(TInt batchsize, TInt ndata){
			this->batchsize = batchsize;
			if (ndata % batchsize == 0){
				this->nsubrounds = ndata/batchsize;
			}
			
			else{
				this->nsubrounds = 1 + ndata/batchsize;
			}
			this->lastbatchsize = ndata - batchsize*(this->nsubrounds - 1);
			this->subround = 0;
			
			
			this->initialising_batch_size = std::min(this->batchsize, ndata);
			this->nchanges_on_batch = std::vector<TInt> (this->nsubrounds, 0);			
		}
};

//template <typename TInt, typename TFloat>
//void set_summaries_minibatch(cluster::BaseCluster<TInt, TFloat> & basecluster, const minibatchapp::MiniBatchApp<TInt> & mba){
	
			



		
}


	
	
#endif
		
