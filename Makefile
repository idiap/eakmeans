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

#################################
#                               #
# USER : PLEASE CONFIGURE HERE  #
#                               # 
#################################

#would you like to use blas? YES/NO
export WITHBLAS=NO

#if you answered yes to the above, 
#the path to the blas library 
#(for example where libopenblas.so is) 
#if not on search path already
export LIBBLASDIR=
#/....../openblas/lib

#where is the header (cblas.h) ?
#if not on searchpath by default
export INCBLASDIR=
#/....../openblas/include

# where to put all .so files other than
#the .so file for python, (wherever you like)  
export LIBDIR=$(shell pwd)

# where to put kmeans.so, the python library
# (wherever you like)
export PYLIBDIR=$(shell pwd)


#################
#   IMPORTANT   #
#################

#you need to 
#add LIBDIR to LD_LIBRARY_PATH
#and 
#add PYLIBDIR to PYTHONPATH
#before building. (ie now)
#for example on linux or mac, 
#in your .bashrc file, do the following now:
#export LD_LIBRARY_PATH=/absolute/path/of/LIBDIR:${LD_LIBRARY_PATH}
#export PYTHONPATH=/absolute/path/of/PYLIBDIR:${PYTHONPATH}

#############################################
#                                           #
# THANKS USER, YOU CAN IGNORE WHAT FOLLOWS  #
#                                           # 
#############################################




all: libraries bins

export CXX :=   g++
export CXXFLAGS := -std=c++11 -O3  -Wall -pedantic -fPIC
export LINKER := g++ #-Wl,--no-as-needed



#export CXX :=   clang++
#export CXXFLAGS := -std=c++11 -O3  -Wall -pedantic -fPIC
#export LINKER := clang++ -std=c++11



#we use this flag to tell the called makefiles that they have been asked to move and shake from from elsewhere
export CALLFROMABOVE=YES

#the directory in which this makefile resides
export ROOTLIBDIR= $(shell pwd)





CURRENTDIR=$(shell pwd)
export UTILITYDIR=${ROOTLIBDIR}/utility
export BINSDIR=${ROOTLIBDIR}/bins
export CYDIR=${ROOTLIBDIR}/cython

LIBRARYDIRS=
BINNAMES=
BINDIRS=

#add blasutil to lib list if USEBLAS is YES, how to make lib and bin/test
export BLASUTILDIR=${UTILITYDIR}/blasutil
ifeq (${WITHBLAS}, YES)
	LIBRARYDIRS+= ${BLASUTILDIR} 
endif
blasutil: 
	@(cd ${BLASUTILDIR} && $(MAKE) all || (echo error making blasutil && exit 41));

#add optionsutil to lib list, how to make lib and bin/test
export OPTIONSUTILDIR=${UTILITYDIR}/optionsutil
LIBRARYDIRS+= ${OPTIONSUTILDIR}
optionsutil: 
	@(cd ${OPTIONSUTILDIR} && $(MAKE) all || (echo error making optionsutil && exit 41));

#add stdthreadutil to lib list, how to make lib and bin/test
export STDTHREADUTILDIR=${UTILITYDIR}/stdthreadutil
LIBRARYDIRS+= ${STDTHREADUTILDIR}
stdthreadutil: 
	@(cd ${STDTHREADUTILDIR} && $(MAKE) all || (echo error making stdthreadutil && exit 41));
	

#add stringutil to lib list, how to make lib and bin/test
export STRINGUTILDIR=${UTILITYDIR}/stringutil
LIBRARYDIRS+=${STRINGUTILDIR}
stringutil: 
	@(cd ${STRINGUTILDIR} && $(MAKE) all || (echo error making stringutil && exit 41));	

#make bin/test for arrutilv2
export ARRUTILV2DIR=${UTILITYDIR}/arrutilv2
arrutilv2: 
	@(cd ${ARRUTILV2DIR}/test && $(MAKE) all || (echo error making arrutilv2 && exit 41));

#make bin/test for randomutil	
export RANDOMUTILDIR=${UTILITYDIR}/randomutil
randomutil: 
	@(cd ${RANDOMUTILDIR}/tests && $(MAKE) all || (echo error making randomutil && exit 41));
	
#make bin/test for sortutil	
export SORTUTILDIR=${UTILITYDIR}/sortutil
sortutil: 
	@(cd ${SORTUTILDIR}/tests && $(MAKE) all || (echo error making randomutil && exit 41));

#add pllkmeans to lib list, how to make lib and bin/test
export PLLKMEANSLDIR=${UTILITYDIR}/pllkmeans
LIBRARYDIRS+=${PLLKMEANSLDIR}
pllkmeans: 
	@(cd ${PLLKMEANSLDIR} && $(MAKE) all || (echo error making pllkmeans && exit 41));

#add kmeans to lib list, how to make lib and bin/test
export KMEANSBDIR=${BINSDIR}/kmeans
BINDIRS+=${KMEANSBDIR}
kmeans: pllkmeans
	@(cd ${KMEANSBDIR} && $(MAKE) all || (echo error making kmeans && exit 41));


export CYKMEANSDIR=${CYDIR}/kmeans
pythonkmeans:
	@(cd ${CYKMEANSDIR} && $(MAKE) all || (echo error making python kmeans module && exit 41));

#	make libraries. will run all makefiles in x \in LIBRARYSOURCEDIRS
# basically, for x \in LIBRARYSOURCEDIRS : @(cd x && $(MAKE))
# (try to change into x and make there, else exit with an error)
libraries: ${LIBRARYDIRS}

	echo ${LIBRARYDIRS}
	@echo making libraries; \
	for dir in ${LIBRARYDIRS}; do \
		echo making lib $$dir; \
		(cd $$dir && $(MAKE) lib); \
	done


# basically, for x \in LIBRARYBINDIRS : @(cd x && $(MAKE))
bins: ${BINDIRS} kmeans
	echo ${BINDIRS}
	@echo making bins; \
	for dir in ${BINDIRS}; do \
		echo making bin $$dir; \
		(cd $$dir && $(MAKE) all); \
	done



# 

SOURCEDIRS=${CYKMEANSDIR} ${LIBRARYDIRS} ${BINDIRS}
clean:
	echo ${SOURCEDIRS}
	@echo making clean; \
	for dir in ${SOURCEDIRS}; do \
		echo cleaning $$dir; \
		(cd $$dir && $(MAKE) clean); \
	done


remove:
	echo ${SOURCEDIRS}
	@echo making clean; \
	for dir in ${SOURCEDIRS}; do \
		echo cleaning $$dir; \
		(cd $$dir && $(MAKE) remove); \
	done


