#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

##########################################################
#compile with blas?
USEBLAS = YES
export LIBBLASDIR=/idiap/user/jnewling/openblas/lib
export INCBLASDIR=/idiap/user/jnewling/openblas/include
##########################################################


CXX :=  g++
CXXFLAGS := -std=c++11 -O3  -Wall -pedantic -fPIC
LINKER := g++ #-Wl,--no-as-needed
LFLAGS := -lpthread
ifeq ($(USEBLAS), YES)
	LFLAGS := ${LFLAGS} -lopenblas -L${LIBBLASDIR}
	CXXFLAGS := ${CXXFLAGS} -D WITHBLAS
	INCLUDEPATHS = -I${INCBLASDIR}
	TARGET := withblaskmeans
	LIBNAME := libwithblaskmeans
	export WITHBLAS
else
	TARGET := blaslesskmeans
	LIBNAME := libblaslesskmeans
endif


SOURCES  := $(wildcard src/*.cpp)
HEADERS	 := $(wildcard src/*.h)
OBJECTS  := $(SOURCES:src/%.cpp=obj/%.o)
OBJECTS_FORLIB := $(filter-out obj/main.o,$(OBJECTS))

all : main lib pythonkmeans


main :  $(OBJECTS)
	@mkdir -p bin	
	$(LINKER) -o bin/${TARGET} $(OBJECTS) $(LFLAGS)
	@echo "Linking for main of ${NAME} done!"

lib : $(OBJECTS_FORLIB)
	@mkdir -p lib
	$(CXX) -shared $^ -o lib/$(LIBNAME).so
	@echo "Shared library ${LIBNAME} made!"

pythonkmeans:
	python setup.py build_ext -b lib

obj/%.o : src/%.cpp  $(HEADERS)
	@mkdir -p obj
	$(CXX) -c  $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@
	@echo "compiled "$<" successfully!"


.PHONEY: clean


clean:
	rm -f  $(OBJECTS)
	@echo "cleanup done!"
	
.PHONEY: remove

remove: clean
	-rm cythonsrc/kmeans.cpp
	-rm -rf build
	-rm -rf bin
	-rm -rf lib
	-rm -rf obj
	@echo "should be 100% clean!"

