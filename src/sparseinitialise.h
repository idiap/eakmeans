
#ifndef SPARSEINITIALISE_H
#define SPARSEINITIALISE_H


#include <memory>
#include <vector>
namespace kmeans{
namespace sparseinit{

//get copyindices guaranteeing that all distinct.
template <typename TFloat, typename TInt>
std::tuple<std::unique_ptr<TFloat []>, std::vector<TInt> > get_initialisation_indices(TInt ncentroids, const sparse::SparseData<TInt, TFloat> & data, TInt data0 = 0, TInt data1 = 0){
	
	if (data1 == 0){
		data1 = data.ndata;
	}


	
	
	TInt ndata_range = data1 - data0;
	std::vector<TInt> initialisation_indices (ncentroids);
	std::unique_ptr<TFloat []> C(new TFloat [ncentroids*data.dimension]);
	
	TInt nattempts = 0;
	TInt currentindex = 0;
	

			
	while (currentindex < ncentroids && nattempts < 5*ncentroids){
		TInt proposal = data0 + rand() % ndata_range;
		
		bool rejected = false;
		for (TInt ci = 0; ci < currentindex; ++ci){
			TFloat l22_diff = sparse::get_l22(

				data.starts[initialisation_indices[ci] + 1] - data.starts[initialisation_indices[ci]],
				data.indices.data() + data.starts[initialisation_indices[ci]],
				data.values.data() + data.starts[initialisation_indices[ci]],

				data.starts[proposal + 1] - data.starts[proposal],
				data.indices.data() + data.starts[proposal],
				data.values.data() + data.starts[proposal]
				
			);
			
			if (l22_diff < 1e-5){
				rejected = true;
				break;
			}
		}
		++nattempts;
		if (rejected == false){
			sparse::todense::zero_and_copy(proposal, data, C.get() + currentindex*data.dimension);
			initialisation_indices[currentindex] = proposal;
			++currentindex;
		}
		else{
		}
	}

	if (currentindex != ncentroids){
		throw std::runtime_error("Tried to find a set of distinct datapoints, but failed (nattempts/ncentroids = 5)"); 
	}

	
	return std::make_tuple (std::move(C), std::move(initialisation_indices));
}



template <typename TFloat, typename TInt>
std::tuple<std::unique_ptr<TFloat []>, std::unique_ptr<TFloat []>, std::unique_ptr<TInt []>, TFloat > 
get_kmeanspp_initialisation(TInt ncentroids, const sparse::SparseData<TInt, TFloat> & data, TInt ind0, TInt ind1){
	throw std::runtime_error("sparse kmeans ++ not yet implemented. Look for inspiration in dense version. Probably common code to be extracted."); 
	
	
	return std::make_tuple (std::unique_ptr<TFloat []> {}, 
	std::unique_ptr<TFloat []> {},
	std::unique_ptr<TInt []> {},
	TFloat {});

}


}
}
#endif

