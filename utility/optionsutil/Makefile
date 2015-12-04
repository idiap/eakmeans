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

NAME:=optionsutil
LIBNAME:=lib${NAME}

INCLUDEPATHS := -I./include
LFLAGS :=

SRCDIR   := src
BINDIR   := bin
INCDIR   := include
OBJDIR   := obj

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
HEADERS	 := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# the library should contain all .o files other than main.o
OBJECTS_FORLIB := $(filter-out $(OBJDIR)/main.o,$(OBJECTS))

all : checkcall main lib

checkcall:
ifneq (${CALLFROMABOVE}, YES)
	@echo "currently there is no support for directly calling this makefile, it should be called by parent makefile (or else, LIBNAME and other variables set in parent should be set here)"
	@exit 43
endif
	
#To make the main executable, link all the object files 
main :  $(OBJECTS)
	$(LINKER) -o bin/$@ $(OBJECTS) $(LFLAGS)
	@echo "Linking for main of ${NAME} done!"

lib : $(OBJECTS_FORLIB)
	$(CXX) -shared $^ -o $(LIBDIR)/$(LIBNAME).so
	@echo "Shared library ${LIBNAME} made!"


$(OBJDIR)/%.o : $(SRCDIR)/%.cpp  $(HEADERS)
	$(CXX) -c  $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@
	@echo "compiled "$<" successfully!"

.PHONEY: clean

clean:
	rm -f  $(OBJECTS) 
	rm -f  $(LIBDIR)/$(LIBNAME).so
	@echo "cleanup done!"

.PHONEY: remove

remove: clean
	rm -f  $(BINDIR)/main
	@echo "executable removed!"

