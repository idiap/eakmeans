from distutils.core import Extension, setup
from Cython.Build import cythonize

import os

#X_library_dir = "./"
openblaslibdir = os.environ["LIBBLASDIR"] #/idiap/user/jnewling/openblas/lib"
	
#LIBBLASDIR	

libname = "kmeans"
	
if "WITHBLAS" not in  os.environ.keys():
	librariestouse = ["blaslesskmeans"]
	print "will build the python library pkmeans without blas (building with blas will make it faster)"

else:
	librariestouse = ["withblaskmeans", "openblas"]
	print "will build the python library pkmeans with blas (good choice)"

#TODO : sort out libname vs "kmeans" below

########## Using Cython directly ###################################
ext = Extension(libname, sources = [os.path.abspath("cythonsrc/kmeans.pyx")], include_dirs = ["cythonsrc", "src"], library_dirs = [openblaslibdir, "lib"], libraries = librariestouse, language = "c++")
setup(name = libname, ext_modules = cythonize(ext))








############ Using precompiled cpp file (no need for Cython) ##########
#ext = Extension("kmeans", sources = [os.path.abspath("./precythonised/kmeans.cpp")], include_dirs = [X_include_dir], library_dirs = [X_library_dir, openblaslibdir], libraries = librariestouse, language = "c++")
#setup(name = "kmeans", ext_modules = [ext])





#if hostname == "goudurix12":
	#X_include_dir = "/idiap/home/jnewling/libraries/utility/pllkmeans/include"
	#X_library_dir = "/idiap/user/jnewling/own/templib"
	#openblaslibdir = "/idiap/user/jnewling/openblas/lib"

#else:
	#X_include_dir = "/home/james/libraries/utility/pllkmeans/include"
	#X_library_dir = "/home/james/oak/own/templib"
	#openblaslibdir = "/home/james/oak/openblas/lib"
