# 
# Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <james.newling@gmail.com>
# All rights reserved.
# 
# eakmeans is a library for exact and approximate k-means written in C++ and
# Python. This file is part of eakmeans. See file COPYING for more details.
# 
# This file is part of eakmeans.
# 
# eakmeans is free software: you can redistribute it and/or modify
# it under the terms of the 3-Clause BSD Licence. See
# https://opensource.org/licenses/BSD-3-Clause for more details.
# 
# eakmeans is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
# COPYING for more details.
# 
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
