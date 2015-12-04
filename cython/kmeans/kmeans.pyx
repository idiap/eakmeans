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

import numpy as np
from libcpp.string cimport string
from libcpp cimport bool
cimport cython
cimport cython.floating
	
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


	return {"C": np.array(C), "L": np.array(L), 'I': np.array(inds0), 'duration':duration, 'niterations':niterations, 'mse':mse, 'output': verbosestring }



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
	
	
	

  	
def get_clustering(X, n_clusters, max_iter = 1e10, n_init = 1, init = "kmeans++", max_time = 1e20, n_threads = 1, seed = None, verbose = 1, algorithm = 'auto', minibatchsize = 100, X_validation = None, validation_period = 5, capture_verbose = False):
	r"""
Input
-------------------
X
	numpy.ndarray of dimension 2 and dtype \in {numpy.float32, numpy.float64}
	the input data. rows considered as datapoints, columns considered as features

n_clusters
	integer type
	the "k" in k-means

algorithm
	string
	the algorithm to use, one of 
	{auto, simple, 
	ham, ann, expSN, expNS, 
	syinSN, syinNS, yin, 
	selkSN, selkNS, elkSN, elkNS}
	where auto will select:
	expNS		if dimension < 6
	syinNS	if 6 <= dimension 40
	selkNS	if 40 <= dimension
	see (TODO my paper) for details
	
	
max_iter 
	integer type
	maximum allowed iterations before returning
	TODO : is mse returned reliable if return due to max_iter?
	
init
	string OR 2-d numpy float array OR 1-d integer array
	if "uniform": centroids initialised as using uniform sampling from X
	if "k-means++": centroids initialised using sampling described in Arthur, D. and Vassilvitskii, S. (2007)
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
	
minibatchsize (UNPUBLISHED)
	integer type
	if algorithm is minibatch, then this sets the batch size

verbose:
	integer type
	if 0: silent
	if 1: print initial and final statistics and number of label changes in each round
	if 2: print statistics every round

capture_verbose:
	bool
	if True : output stream is put into `output'
	else : output goes directly to terminal
	
X_validation:
	TODO
	
validation_period:
	TODO

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
	
	
	
	#TODO : early error catching here


	
	ncentroids = n_clusters
	ndata, dimension = X.shape
	float_type = X.dtype
	cout_verbosity = verbose

	if algorithm == 'auto': 
		if dimension < 6:
			algorithm = 'expNS'
		elif dimension < 40:
			algorithm = 'syinNS'
		else:
			algorithm = 'selkNS'
	
	
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
			data_indices_init_from = init
	
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
		
	return dangerwrap(lambda : basekmeans(algorithm, n_threads, ndata, dimension, X.ravel(), ncentroids, cout_verbosity, initialisation_method, C_init.ravel(), data_indices_init_from.ravel(), setseed, seed, max_iter, max_time, minibatchsize, nvaldata, X_validation.ravel(), validation_period, capture_verbose))
