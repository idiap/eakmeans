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

"""
Q: What is the effect on runtime of verbose output?
"""

import sys
sys.path.append(os.path.abspath(".."))
sys.path.append(os.path.abspath("../.."))


import pyximport 
pyximport.install(reload_support = True)

import kmeans
reload(kmeans)

import numpy as np
import numpy.random as npr

import time
import sys
	
from IPython.core.debugger import Tracer



def on_syinNS():
	datafile = None # path to data
	data = np.loadtxt(datafile)
	ndata, dim = data.shape
	data_validation = None
	

	for ncentroids in [10, 100, 1000]:
		C0 = npr.permutation(data)[0:ncentroids, :]
		print "with #centroids = ", ncentroids
		print "verbosity \ttime"
		for verb, alg in zip([0,1,2], ['syinSN', 'syinSN', 'syinSN']):
			t0 = time.time()	
			X = kmeans.get_clustering(X = data, n_clusters = ncentroids, init = C0, algorithm = alg, n_threads = 1, verbose = verb, X_validation = data_validation, validation_period = 1, capture_verbose = True)
			print verb, "\t\t", time.time() - t0
			

def on_standard():
	data = npr.randn(20000, 5)
	ndata, dim = data.shape
	data_validation = None
	
	for ncentroids in [10, 100, 1000]:
		C0 = npr.permutation(data)[0:ncentroids, :]
		print "with #centroids = ", ncentroids
		print "verbosity \ttime"
		alg = 'simple'
		for verb, alg in zip([0,1,2], [alg, alg, alg]):
			t0 = time.time()	
			X = kmeans.get_clustering(X = data, n_clusters = ncentroids, init = C0, algorithm = alg, n_threads = 1, verbose = verb, X_validation = data_validation, validation_period = 1, capture_verbose = True)
			print verb, "\t\t", time.time() - t0
		
		
		
		
	
