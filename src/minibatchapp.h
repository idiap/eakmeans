/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


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
		
