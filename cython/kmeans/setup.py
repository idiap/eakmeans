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

from distutils.core import Extension, setup
from Cython.Build import cythonize

import os

if 'CALLFROMABOVE' not in os.environ.keys():
	print "CALLFROMABOVE NOT SET, EXITING..."

else:
	X_include_dir = os.environ['PLLKMEANSLDIR'] + "/include"
	X_library_dir = os.environ['LIBDIR'] 
	openblaslibdir = os.environ['LIBBLASDIR']



librariestouse = ["pllkmeans", "stringutil", "stdthreadutil"]

if os.environ["WITHBLAS"] == "YES":
	librariestouse += ["openblas", "blasutil"]

########### Using Cython directly ###################################
#ext = Extension("kmeans", sources = [os.path.abspath("kmeans.pyx")], include_dirs = [X_include_dir], library_dirs = [X_library_dir, openblaslibdir], libraries = librariestouse, language = "c++")
#setup(name = "kmeans", ext_modules = cythonize(ext))


########### Using precompiled cpp file (no need for Cython) ##########
ext = Extension("kmeans", sources = [os.path.abspath("./precythonised/kmeans.cpp")], include_dirs = [X_include_dir], library_dirs = [X_library_dir, openblaslibdir], libraries = librariestouse, language = "c++")
setup(name = "kmeans", ext_modules = [ext])



