import kmeans

import numpy as np
import numpy.random as npr

old_seed = npr.randint(100000)

ndata = 10000
dimension = 300
n_clusters  = 10
npr.seed(1012)

X = npr.randn(ndata, dimension)
C0 = 1.001*npr.randn(n_clusters, dimension)

npr.seed(old_seed)

bla = kmeans.get_clustering(X = X, n_clusters = n_clusters, init = "BF", verbose = 1, seed = old_seed)
