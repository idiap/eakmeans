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

Making and loading all things MNIST (without labels) for clustering tasks.
See pythongold if you're interested in labeled data.

* Standard mnist dataset
* Enlarged mnist datasets
* 00 -> 99 mnist datasets (concatenated)
* Random projections of all above (projected) (TODO for concatenated still)

"""

import os, struct
from array import array
from cvxopt.base import matrix
import numpy as np
import socket
import commands

from IPython.core.debugger import Tracer
import numpy as np
import numpy.random as npr

#set the path to infiniminst
infipath = '/idiap/temp/jnewling/data/infinitemnist/infimnist'
infiexec = os.path.join(infipath, "infimnist")

#set path to data dirs

paths = {
'original' : os.path.join(infipath, 'generateddata/original'), 
'projected' : os.path.join(infipath, 'generateddata/projected'),
'concatenated' : os.path.join(infipath, 'generateddata/concatenated')}

if (socket.gethostname() == 'tamarind'):
	cooki = "/home/james/data/mnist"


def make_MNIST(dataset = "original", ndata = 100000, dimension = None):

  print "in make MNIST"
  if dataset == "original":
    cwd = os.getcwd()
    os.chdir(infipath)
    datafn = os.path.join(paths['original'], 'original-%d-ubyte'%(ndata,))
    print commands.getstatusoutput('infimnist pat 0 %d > %s'%(ndata, datafn))
    os.chdir(cwd)
    
  if dataset == "projected":
    datafn = os.path.join(paths['projected'], 'projected_%d_%d.npy'%(ndata,dimension))
    datafn_preprojected = os.path.join(paths['original'], 'original-%d-ubyte'%(ndata,))
    if not os.path.exists(datafn_preprojected):
      make_MNIST(dataset = "original", ndata = ndata)
      X = read_MNIST("original", ndata)
      os.remove(datafn_preprojected)
    else:
      print "reading original"
      X = read_MNIST("original", ndata)
    
    npr.seed(1011) 
    projection_matrix = npr.randn(784, dimension)
    X_proj = np.dot(X, projection_matrix).astype(np.float32)
    np.save(arr = X_proj, file=  datafn)      
        
def get_two_digs(X, ndata):
  n_orig, dimension = X.shape
  twod = np.empty((ndata, 2*dimension), dtype = np.float32)
  twod_left = twod[:, 0:dimension]
  twod_right = twod[:, dimension:2*dimension]
  for i in range(n_orig):
    if i*n_orig < ndata:
      twod_left[i*n_orig:min((i + 1)*n_orig, ndata)] = X[i].reshape(1, -1)
      nright = min((i + 1)*n_orig, ndata) - i*n_orig
      twod_right[i*n_orig:i*n_orig + nright, :] = X[0:nright, :]

  #Tracer()()
  return twod
 
  
def read_MNIST(dataset = "original", ndata = 100000, dimension = None):
    """
    Python function for importing the MNIST data set.
    """

    if dataset == "original":
      fname_img = os.path.join(paths['original'], 'original-%d-ubyte'%(ndata,))
      if not os.path.exists(fname_img):
        make_MNIST("original", ndata)
      
      
      fimg = open(fname_img, 'rb')
      magic_nr, size, rows, cols = struct.unpack(">IIII", fimg.read(16))
      img = array("B", fimg.read())
      fimg.close()
  
      images =  matrix(0, (ndata, rows*cols))
      for i in range(ndata):
        images[i, :] = img[ i*rows*cols : (i+1)*rows*cols ]
    
      return np.array(images, dtype = np.float32)
      
    elif dataset == "projected":
      datafn = os.path.join(paths['projected'], 'projected_%d_%d.npy'%(ndata,dimension))
      if not os.path.exists(datafn):
        make_MNIST("projected", ndata, dimension)
      
      return np.load(datafn)
    
    #concatenated  (never actually saved as such)
    else:
      source_ndata = int(np.sqrt(ndata)) + 1
      sourcefn = os.path.join(paths['original'], 'original-%d-ubyte'%(source_ndata,))
      if not os.path.exists(sourcefn):
        make_MNIST("original", source_ndata)
      
      X = read_MNIST(dataset = "original", ndata = source_ndata)
      return get_two_digs(X, ndata)
      
      
      
    
