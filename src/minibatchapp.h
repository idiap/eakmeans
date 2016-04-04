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
		
