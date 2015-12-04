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

INCLUDEPATHS := -I../include -I${INCBLASDIR}
LFLAGS   := -lopenblas -L${LIBBLASDIR}

SRCDIR   := src
BINDIR   := bin
INCDIR   := ../include
OBJDIR   := obj
EXE := main

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
HEADERS	 := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all : checkcall ${EXE}

checkcall:
ifneq (${CALLFROMABOVE}, YES)
	@echo "currently there is no support for directly calling this makefile, it should be called by parent makefile"
	@exit 43
endif

#To make the main executable, link all the object files 
${EXE} :  $(OBJECTS)
	$(LINKER) -o bin/$@ $(OBJECTS) $(LFLAGS)
	@echo "Linking for ${EXE} of ${NAME} done!"

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp  $(HEADERS)
	$(CXX) -c  $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@
	@echo "compiled "$<" successfully!"

.PHONEY: clean

clean:
	rm -f  $(OBJECTS) 
	rm -f $(LIBDIR)/$(LIBNAME).so
	@echo "cleanup done!"

.PHONEY: remove

remove: clean
	rm -f  $(BINDIR)/$(TARGET)
	@echo "executable removed!"








#SRCDIR   := src
#OBJDIR   := obj
#BINDIR   := bin
#INCDIR   := include
#LIBDIR   := lib

#SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
#HEADERS	 := $(wildcard $(INCDIR)/*.h)
#OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)



##To make the main executable, link all the object files 
#main: $(OBJECTS)
	#@$(LINKER) $(LFLAGS) $(OBJECTS) -o ${BINDIR}/$@
	#@echo "Linking done!"

#$(OBJDIR)/%.o : $(SRCDIR)/%.cpp  $(HEADERS)
	#@$(CXX) -c $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@
	#@echo "compiled "$<" successfully!"

#.PHONEY: clean

#clean:
	#@rm -f $(OBJECTS)
	#@echo "cleanup done!"

#.PHONEY: remove

#remove: clean
	#@rm -f $(BINDIR)/main
	#@echo "executable removed!"
