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
import numpy as np
import numpy.random as npr
import time
from IPython.core.debugger import Tracer

import sys
sys.path.insert(0, "../lib")

import kmeans
reload(kmeans)


def example_1(ndata = 1e4, dimension = 100, dtype = np.float64):
	"""
	basic use : cluster random data
	"""
	data = npr.randn(ndata, dimension).astype(dtype)
	clustering = kmeans.get_clustering(X = data, n_clusters = 60, algorithm = 'auto', verbose = 2)
	

def example_2():
	"""
	compare algorithms which may be good in low-d 
	ham, ann, expSN, expNS, syinSN, syinNS, yin,
	on dataset ldfpads.txt (~160000 points in 3 dimensions)
	"""
	data = np.loadtxt('ldfpads.txt')
	print "Data shape : ", data.shape
	seed = npr.randint(100000)
	algs = ['ham','ann', 'exp-sn', 'exp-ns', 'syin-sn', 'syin-ns', 'yin']

	times = dict.fromkeys(algs)
	for alg in algs:
		clustering = kmeans.get_clustering(X = data, n_clusters = 1000, algorithm = alg, verbose = 1, n_threads = 1, seed = seed)
		times[alg] = clustering['duration']
		
	
	print "TIMES:"
	for alg in algs:
		print alg, " : ", times[alg]
		

def example_3():
	"""
	compare selkNS to standard algorithm on random data. selkNS is faster, but not by as much as when the data has structure.
	"""	
	algs = ['sta','selk-ns']
	times = dict.fromkeys(algs)
	data = npr.randn(50000, 25).astype(np.float64)
	seed = 1011
	for alg in algs:
		clustering = kmeans.get_clustering(X = data, n_clusters = 1000, algorithm = alg, verbose = 1, n_threads = 1, seed = seed)
		times[alg] = clustering['duration']
	

	print "TIMES:"
	for alg in algs:
		print alg, " : ", times[alg]
		

def example_4():
	"""
	compare to scikitlearn implementation of kmeans
	"""

	import sklearn.cluster as skc
	import time
	
	ndata = 50000
	dimension = 10
	ncentroids = 1000
	data = npr.randn(ndata, dimension).astype(np.float64)

	centroids0 = data[0:ncentroids, :]

	t0 = time.time()
	kmeans.get_clustering(X = data, init = centroids0, n_clusters = ncentroids, algorithm = 'auto', verbose = 1, n_threads = 1)
	t1 = time.time()

	sklearner = skc.k_means(X = data, n_clusters = ncentroids, max_iter = 1000, n_init = 1, init = centroids0, precompute_distances = False, verbose = True, n_jobs = 1, return_n_iter = True, tol = 0.0)
	t2 = time.time()	
	
	kmeans_time = t1 - t0
	sklearner_time = t2 - t1
	
	print "sklearn : ", sklearner_time, " s"
	print "this kmeans: ",  kmeans_time, " s"
			
	

