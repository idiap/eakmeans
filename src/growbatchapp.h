#ifndef GBAPP_H
#define GBAPP_H

#include <memory>

namespace growbatchapp{
	
template <typename TInt, typename TFloat>
class GBApp {
	public:
	
	
		/* amount of data which is active. Initially determined by user, thereafter grows by factor of growthfactor when nec */
		TInt ndata_active; 

		/* Hacky variable for printing purposes */
		TFloat d_C__over__d_AB;
		
		/* will be 2 */
		TFloat growthfactor; 
		
		/* definition depends on class, for Grow Batch Partitional: if \|C_A - C_B\|_2 > threshold * \|C_{t} - C_{t-1}\|, then grow by growthfactor. will be 1.0 */
		TFloat threshold; 
		
		/* amount of data which was active in previous round. Either ndata_active or ndata_active/2. */
		TInt ndata_active_previous; 

		/* used to determine if exapansion should take place (and maybe other things) */
		std::unique_ptr<TFloat []> delta_C;


};

template <typename TInt, typename TFloat>
class GBMseApp {
	public:
		std::vector<TFloat> sse_by_cluster;
		std::vector<TFloat> mse_by_cluster;
		std::unique_ptr<TFloat []> dn; //distance to nearest. TODO: checl that not this quantity squared.	
};

}






#endif
		
