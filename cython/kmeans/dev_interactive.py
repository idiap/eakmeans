#EAKMeans is a fast Exact K-means library written in C++ with 
#command-line interface, shared library + header files and 
#Python bindings

#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#This file is part of EAKMeans.

#EAKMeans is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as
#published by the Free Software Foundation.

#EAKMeans is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.

#for interactive loading (a change in .pyx is detected and cython is rerun)
import pyximport 
pyximport.install(reload_support = True)

import kmeans
reload(kmeans)

import sys
sys.path.append("../../python/datasets/")
import numpy as np
import numpy.random as npr


#import mnist 
#reload(mnist)

import time
import sys
	
from IPython.core.debugger import Tracer
	



nbatches = 1
ndata = 6000

dimension = 10
data =  npr.randn(ndata, dimension)

#data =  mnist.read_MNIST(dataset = "projected", ndata = ndata, dimension = dimension).astype(np.float64)


ncentroids = 20

C0 = data[0:ncentroids,:]
print "entering clustering"

#TODO : `change' minibatchsize to GBsize0

#reso = kmeans.get_clustering(X = data, n_clusters = ncentroids, init = C0, algorithm = "simple", n_threads = 1, verbose = 2, minibatchsize = ndata/2, X_validation = data, validation_period = 1, capture_verbose = False, max_iter = 5)

reso = kmeans.get_clustering(X = data, n_clusters = ncentroids, init = C0, algorithm = "gbsimple", n_threads = 1, verbose = 2, minibatchsize = ndata/2, X_validation = data, validation_period = 1, capture_verbose = False, max_iter = 15)

