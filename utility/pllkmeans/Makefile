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

USEGPROF:=NO

NAME:=pllkmeans
LIBNAME:=lib${NAME}

#ALL .h files listed here,
INCLUDEPATHS := -I./include -I./include/kmeans -I./include/kmeans/X_step_algorithms -I./include/kmeans/exact -I./include/kmeans/minibatch  -I./include/kmeans/growbatch -I${UTILITYDIR}/optionsutil/include -I${UTILITYDIR}/stdthreadutil/include -I${UTILITYDIR}/randomutil/include  -I${UTILITYDIR}/stringutil/include -I${CLUSTERINGUTILDIR}/include -I${UTILITYDIR}/sortutil/include -I${UTILITYDIR}/blasutil/include -I${UTILITYDIR}/arrutilv2/include

LFLAGS   :=  -lpthread  -L${LIBDIR} -lstdthreadutil -lstringutil  

ifeq ($(WITHBLAS), NO)
	PROJECT:=blaslesskmeans_test
	INCLUDEPATHS :=  -I${UTILITYDIR}/arrutilv2/blaslessbase/include ${INCLUDEPATHS}  
	OBJDIR := obj/blaslessobj
else
	PROJECT:=withblaskmeans_test
	INCLUDEPATHS := ${INCLUDEPATHS} -I${UTILITYDIR}/blasutil/include  -I${UTILITYDIR}/arrutilv2/withblasbase/include -I${INCBLASDIR}
	OBJDIR := obj/withblasobj
	LFLAGS := ${LFLAGS} -lblasutil -lopenblas	-L${LIBBLASDIR}  
endif

ifeq (${USEGPROF}, YES)
	LFLAGS += -pg
endif

SRCDIR   := src
BINDIR   := tests
INCDIR   := include

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
HEADERS	 := $(wildcard $(INCDIR)/*.h)

OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o) 

# the library should contain all .o files other than ${PROJECT}.o
OBJECTS_FORLIB := $(filter-out $(OBJDIR)/${PROJECT}.o, $(OBJECTS)) 

all : checkcall lib ${BINDIR}/${PROJECT}

checkcall:
ifneq (${CALLFROMABOVE}, YES)
	@echo "currently there is no support for directly calling this makefile, it should be called by parent makefile (or else, LIBNAME and other variables set in parent should be set here)"
	@exit 43
endif

#To make the main executable, link all the object files 
${BINDIR}/${PROJECT} : $(OBJECTS)
	$(LINKER) -o $@ $(OBJECTS) $(LFLAGS)
	@echo "Linking for ${BINDIR}/${PROJECT} of ${NAME} done!"

lib : $(OBJECTS_FORLIB)
	$(CXX) -shared $^ -o $(LIBDIR)/$(LIBNAME).so
	@echo "Shared library ${LIBNAME} made!"


#set up some dependencies which should be recognised while prototyping:
#how to make pllkmeansfuncs.o -- note the use of -O3 flag, see http://stackoverflow.com/questions/6726327/does-extern-template-prevent-inlining-of-functions and ask what is going on?
PLLKMEANSFUNCSHEADERS = $(addprefix ${INCDIR}/, pllkmeansfuncs.hpp pllkmeansfuncs_void.h pllkmeansfuncs_nonvoid.h)
$(OBJDIR)/pllkmeansfuncs.o : $(SRCDIR)/pllkmeansfuncs.cpp ${PLLKMEANSFUNCSHEADERS} $(INCDIR)/pllcluster.h
	$(CXX) -c  $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@
	@echo "(SRC) compiled "$<" successfully!"
	

#how to build an object file, generic
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp 
	$(CXX) -c  $(CXXFLAGS)  $(INCLUDEPATHS) $< -o $@
	@echo "(SRC) compiled "$<" successfully!"

.PHONEY: clean

clean:
	rm -f $(wildcard $(OBJDIR)/*.o)	
	rm -f $(LIBDIR)/$(LIBNAME).so
	@echo "cleanup done!"

.PHONEY: remove

remove: clean
	rm -f  $(BINDIR)/main
	@echo "executable removed!"

