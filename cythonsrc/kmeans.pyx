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
from libcpp.string cimport string
from libcpp cimport bool
cimport cython
cimport cython.floating
import random
	
cdef extern from "pllkmeansfuncs_void.h" namespace "cluster":	
	void v_solveiolessf(const string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const float * const data, size_t ncentroids, int cout_verbosity, const string & initialisation_method, const float * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, float maxtime, size_t maxrounds, float * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t &  niterations, float & mse, size_t minibatchsize, size_t nvaldata, const float * const valdata, size_t valperiod, bool captureverbose, string & verbosestring) except + #except + as it might throw an exception, "so catch it cython"

	
	void v_solveiolessd(const string & algorithm, size_t nthreads, size_t ndata, size_t dimension, const double * const data, size_t ncentroids, int cout_verbosity, const string & initialisation_method, const double * const C_init, const size_t * const data_indices_init_from, bool setseed, size_t seed, double maxtime, size_t maxrounds, double * const C, size_t * const L, size_t * const inds0, size_t & duration, size_t & niterations, double & mse, size_t minibatchsize, size_t nvaldata, const double * const valdata, size_t valperiod, bool captureverbose, string & verbosestring) except + 
	

		
def basekmeans(algorithm, nthreads, ndata, dimension, cython.floating [:] X, ncentroids, cout_verbosity, initialisation_method, cython.floating [:] C_init, size_t [:] data_indices_init_from, setseed, seed, max_iter, max_time, minibatchsize, nvaldata, cython.floating [:] X_validation, validation_period, bool capture_verbose):



	cdef size_t [:] L = np.empty((ndata,), dtype = np.uint64)
	cdef size_t [:] inds0 = np.empty((ncentroids,), dtype = np.uint64)
	cdef size_t duration = 0
	cdef size_t niterations = 0
	cdef cython.floating mse = 0
	
	cdef void (*grail)(int, char*) except +
		
	cdef void (*graily)(const string &, size_t, size_t, size_t, const cython.floating * const , size_t, int, const string & , const cython.floating * const, const size_t * const, bool, size_t, cython.floating, size_t, cython.floating * const, size_t * const, size_t * const , size_t &, size_t &, cython.floating &, size_t, size_t, const cython.floating * const, size_t, bool, string &) except +
		
		
	cdef cython.floating [:] C
	if cython.floating is double:
		C = np.empty((ncentroids*dimension,), dtype = np.float64)
		graily = & v_solveiolessd
		
	else:
		C = np.empty((ncentroids*dimension,), dtype = np.float32)
		graily = & v_solveiolessf

	#is it initialised as a nullptr?
	cdef const cython.floating * C_init_ptr = NULL 
	cdef const size_t * data_indices_init_from_ptr = NULL
	cdef const cython.floating * X_validation_ptr = NULL
	
	if initialisation_method == "from_C":
		C_init_ptr = &C_init[0]

	elif initialisation_method == "from_indices":
		data_indices_init_from_ptr = &data_indices_init_from[0]

	if nvaldata > 0:
		X_validation_ptr = &X_validation[0]
	
	cdef string verbosestring = ""	
	
	graily(algorithm, nthreads, ndata, dimension, &X[0], ncentroids, cout_verbosity, initialisation_method,  C_init_ptr, data_indices_init_from_ptr, setseed, seed, max_time, max_iter, &C[0], &L[0], &inds0[0], duration, niterations, mse, minibatchsize, nvaldata, X_validation_ptr, validation_period, capture_verbose, verbosestring)


	return {"C": np.array(C).reshape(ncentroids,dimension), "L": np.array(L), 'I': np.array(inds0), 'duration':duration, 'niterations':niterations, 'mse':mse, 'output': verbosestring }



import multiprocessing as mpr
def dangerwrap(f):
	"""
	I assume f is a function which returns 
	an object and takes no parameters
	"""
	event  = mpr.Event()
	q = mpr.Queue()
	
	def signalling_f():
		try:
			q.put(f())
		except Exception as e:
			print "Caught exception in dangerwrap:"
			print e
			q.put(e)
			event.set()
			return None
			
		event.set()

		
	
	f_process = mpr.Process(target = signalling_f)
	f_process.start()
	try:
		event.wait()
		
	except KeyboardInterrupt:
		f_process.terminate()
		f_process.join()
		raise KeyboardInterrupt("Caught KeyboardInterrupt in dangerwrap")

	return q.get()
	
	
	

		
 # , n_init = 1
 
def get_clustering(X, n_clusters, max_iter = 1e10, init = "kmeans++", max_time = 1e20, n_threads = 1, seed = None, verbose = 1, algorithm = 'auto', minibatchsize = 100, X_validation = None, validation_period = 5, capture_verbose = False):
	"""  
	
Input
-------------------
X
	numpy.ndarray of dimension 2 and dtype of either numpy.float32 or numpy.float64
	the input data. rows considered as datapoints, columns considered as features

n_clusters
	integer type
	the "k" in k-means

algorithm 
	string
	the algorithm to use, the best algorithm for the job is probably one of (see [arXiv:1602.02514] and [arXiv:1602.02934])
	-----------------------------------------------------------------------
	syin-ns : Simplified Yinyang algorithm with ns-bounding, the fastest exact algorithm in dimensions ~ 7->50, memory O(NK) but lighter in memory that selk-ns
	exp-ns : Exponion with ns-bounding, the fastest non-exact algorithm in dimension below ~ 7, memory O(N)
	selk-ns : Simplified version of Elkan with ns-bounding, the fastest exact algorithm in dimensions greater that ~50, memory O(NK)
	turbobatch-rho : Turbobatch algorithm [arXiv:1602.02934], fastest non-exact algorithm, a good choice in general, memory O(NK)
	auto : one of three depending on dimension: ( 1 exp-ns 6 syin-ns 50 selk-ns )	
	------------------------------------------
	other algorithms which are implemented are
	------------------------------------------
	standard : Lloyd's algorithm, the simple exact algorithm
	selk-sn : Simplified version of Elkan with sn-bounding, in general outperformed by selk-ns
	elk-sn : Elkan, in general outperformed by selk-ns
	syin-sn : Simplified Yinyang algorithm with sn-bounding, in general outperformed by syin-ns	
	ham : Hamerly, in general outperformed by exp-ns
	ann : Annulus, in general outperformed by exp-ns
	minibatch-f : Fixed version of minibatch, in general outperformed by turbobatch-rho
	minibatch : Unfixed version of minibatch, in general outperformed by turbobatch-rho
	exp-sn : Exponion with ns-bounding, in general slightly outperformed by exp-ns
	elk-ns : Elkan with ns-bounding, in general outperformed by selk-ns
	yin-sn : Yinyang algorithm, in general outperformed by syin-ns
	growbatch-rho : Growbatch algorithm [arXiv:1602.02934],  in general outperformed by turbobatch-rho O(N)
	
max_iter 
	integer type
	maximum allowed iterations before returning
	TODO : is mse returned reliable if return due to max_iter?
	
init
	string OR 2-d numpy float array OR 1-d integer array
	if "uniform": centroids initialised as using uniform sampling from X
	if "kmeans++": centroids initialised using sampling described in Arthur, D. and Vassilvitskii, S. (2007)
	if "BF": use Bradley and Fayyad (with J = 10, as in Celebi et al.),
	if "BFS": like Bradley and Fayyad, but use initialisation of best partion cluster set, not its final centers
	if 2-d numpy float array : assumed to be initialising centroids
	if 1-d integer array : assumed to be indices of data used to initialise centroids

max_time
	integer type
	maximum allowed in milliseconds time before returning: no new iteration is started after max_time, the most recent results returned

n_threads
	integer type
	number of threads to use
	
seed 
	integer type
	initialising random seed
	
minibatchsize
	integer type
	if algorithm is a variant on minibatch or growbatch, then this sets the (initial) batch size

verbose:
	integer type
	if 0: silent
	if 1: print initial and final statistics and number of label changes in each round
	if 2: print statistics every round

capture_verbose:
	bool
	if True : output stream is put into `output'
	else : output goes directly to terminal
	

Output
-------------------
dict with keys

C
	numpy.ndarray of dimension 2 and dtype the same as X
	the final centroids

L
	numpy.ndarray of dimension 1 and dtype TODO
	final cluster indices of data

I
	numpy.ndarray of dimension 1 and dtype TODO OR None
	if relevant, the indices of data points used to initialise centroids 


duration
	float
	execution time in milliseconds

iterations
	int
	number of iterations

mse
	float
	mean squared distance over data to nearest centroid
	"""
	

	valid_algorithms =	["turbobatch-rho", "exp-ns", "selk-ns", "syin-ns", "elk-sn", "syin-sn", "ham", "ann", "minibatch-f", "standard", "selk-sn","elk-ns", "exp-sn", "yin-sn", "minibatch", "growbatch-rho"]
	
	mapfromusertoraw = {}
	mapfromusertoraw["turbobatch-rho"] =  "gbmse3v1";
	mapfromusertoraw["growbatch-rho"] =  "gbmse";
	mapfromusertoraw["exp-ns"] =  "expNS" ;
	mapfromusertoraw["selk-ns"] =  "selkNS"; 
	mapfromusertoraw["syin-ns"] =   "syinNS";
	mapfromusertoraw["elk-sn"] =   "elkSN";
	mapfromusertoraw["syin-sn"] =   "syinSN";
	mapfromusertoraw["ham"] =  "ham";
	mapfromusertoraw["ann"] =   "ann";
	mapfromusertoraw["minibatch-f"] =  "minibatch";
	mapfromusertoraw["standard"] =  "exactsimplebatch";
	mapfromusertoraw["selk-sn"] =  "selkSN";
	mapfromusertoraw["elk-ns"] =   "elkNS";
	mapfromusertoraw["exp-sn"] =   "expSN";
	mapfromusertoraw["yin-sn"] =  "p17v6";
	mapfromusertoraw["minibatch"] =   "standardminibatch";
	
	ncentroids = n_clusters
	ndata, dimension = X.shape
	float_type = X.dtype
	cout_verbosity = verbose

	if algorithm == 'auto': 
		if dimension < 6:
			algorithm = 'expNS'
		elif dimension < 50:
			algorithm = 'syinNS'
		else:
			algorithm = 'selkNS'
	
	else:
		if algorithm not in valid_algorithms:
			print "algorithm ", algorithm, " is not a valid option. Please choose one of the algorithms suggested"
			return 1
		else:
			algorithm = mapfromusertoraw[algorithm]
	
	if seed == None:
		setseed = False
		seed = 0
		
	else:
		setseed = True

	#A hack: create python objects which can slip through a C++ pointer net 
	data_indices_init_from = np.empty((1,), dtype = np.uint64)
	C_init = np.empty((1,), dtype = float_type)
	
	if isinstance(init, np.ndarray):
		if init.ndim == 2:
			initialisation_method = "from_C"
			C_init = init.ravel()
			
		elif init.ndim == 1:
			initialisation_method = "from_indices"
			data_indices_init_from = init.copy()
			data_indices_init_from.sort()
	
	elif "BF" in init:
		
		random.seed(seed)
		J = 10
		index_shuffle = random.sample(xrange(ndata), ndata)
		index_shuffle_reverse = np.array(np.argsort(index_shuffle), dtype = np.uint64)
    
		X_copy = X[index_shuffle]
		partition_Cs = np.empty((0, dimension), dtype = X.dtype)
		
		#will store indices, accoring to original X:
		partition_true_start_indices = []
		
		for j in range(J):
			partition_start = int(ndata * (j + 0.)/(J + 0.))
			partition_end = int(ndata * (j + 1.)/(J + 0.))
			
			clustering_j = get_clustering(X_copy[partition_start:partition_end, :], n_clusters = n_clusters, max_iter = max_iter, init = "uniform", max_time = max_time, n_threads = n_threads, seed = 1011, verbose = 0, algorithm = 'auto', minibatchsize = minibatchsize, X_validation = None, validation_period = 0, capture_verbose = False)
			indices_in_X_copy = clustering_j['I'] + partition_start
			indices_in_X = index_shuffle_reverse[indices_in_X_copy]
			partition_true_start_indices.append(indices_in_X)
			
			partition_Cs = np.vstack([partition_Cs, clustering_j['C']])
		
#		partition_Cs_list = partition_Cs.tolist()
		potential_initialiser_energies = []
		random.seed(1011)
		for j in range(J):

#			potential_initialisers.append(np.array(random.sample(partition_Cs_list, n_clusters), dtype = X.dtype))

			clustering_j = get_clustering(partition_Cs, n_clusters = n_clusters, max_iter = max_iter, init = np.arange(j*n_clusters, (j+1)*n_clusters, dtype = np.uint64), max_time = max_time, n_threads = n_threads, seed = 1011, verbose = 0, algorithm = 'auto', minibatchsize = minibatchsize, X_validation = X_validation, validation_period = 0, capture_verbose = False)
		
			potential_initialiser_energies.append(clustering_j['mse'])
		
		lowest_energy_index = np.argmin(potential_initialiser_energies)
		

    #indices ?
		if init == "BFS":
			data_indices_init_from = partition_true_start_indices[lowest_energy_index]
			data_indices_init_from.sort()
			initialisation_method = "from_indices"
    

		if init == "BF":
			C_init = partition_Cs[j*n_clusters:(j+1)*n_clusters, :].ravel()
			initialisation_method = "from_C"



		
		
		
		
		
	else:
		initialisation_method = init
	
	if isinstance(X_validation, np.ndarray):
		if X_validation.dtype != X.dtype:
			raise RuntimeError("Type of X_validation and X should be the same")
			
		elif X_validation.ndim != 2:
			runtimeerror = "X_validation should be a 2 dimensional array, it is a %s array."%(X_validation.ndim,)
			raise RuntimeError(runtimeerror)

		nvaldata = X_validation.shape[0]

	else :
		#a hack to let slip through 
		X_validation = np.empty((1,), dtype = float_type)
		nvaldata = 0
	
	if ncentroids < 2:
		raise RuntimeError("number of clusters should be 2 or larger")
	
	if (ncentroids == 2) and ("auto" in algorithm or "exp" in algorithm or "yin" in algorithm):
		raise RuntimeError("Currently, no exp or yin algorithms can be run with ncentroids == 2. Please select another algorithm (not auto/exp/yin) ")
	
	
	return dangerwrap(lambda : basekmeans(algorithm, n_threads, ndata, dimension, X.ravel(), ncentroids, cout_verbosity, initialisation_method, C_init.ravel(), data_indices_init_from.ravel(), setseed, seed, max_iter, max_time, minibatchsize, nvaldata, X_validation.ravel(), validation_period, capture_verbose))

