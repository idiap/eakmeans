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
	$(LINKER) -o bin/${TARGET} $(OBJECTS) $(LFLAGS)
	@echo "Linking for main of ${NAME} done!"

lib : $(OBJECTS_FORLIB)
	$(CXX) -shared $^ -o lib/$(LIBNAME).so
	@echo "Shared library ${LIBNAME} made!"

pythonkmeans:
	python setup.py build_ext -b lib

obj/%.o : src/%.cpp  $(HEADERS)
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
	-rm bin/*
	-rm lib/*
	@echo "should be 100% clean!"

