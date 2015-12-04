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

USEGPROF := NO

NAME:=pllkmeans

INCLUDEPATHS := -I${UTILITYDIR}/optionsutil/include -I${UTILITYDIR}/pllkmeans/include
LFLAGS   :=  -lpthread  -L${LIBDIR} -loptionsutil -lpllkmeans -lstdthreadutil -lstringutil

ifeq ($(WITHBLAS), NO)
	PROJECT:=blaslesskmeans
	OBJDIR   := obj/blaslessobj
else
	PROJECT:=withblaskmeans
	OBJDIR   := obj/withblasobj
	LFLAGS := ${LFLAGS} -lblasutil -lopenblas	-L${LIBBLASDIR}
endif

ifeq (${USEGPROF}, YES)
	LFLAGS += -pg
endif

SRCDIR   := src
BINDIR   := bin

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

all : checkcall main

checkcall:
ifneq (${CALLFROMABOVE}, YES)
	@echo "currently there is no support for directly calling this makefile, it should be called by parent makefile"
	@exit 43
else 
# all is fine
endif	

main : ${BINDIR}/${PROJECT}

${BINDIR}/${PROJECT} :  $(OBJECTS)
	$(LINKER) -o $@ $(OBJECTS) $(LFLAGS)
	@echo "Linking for ${BINDIR}/${PROJECT} of ${NAME} done!"

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp 
	$(CXX) -c  $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@
	@echo "compiled "$<" successfully!"

.PHONEY: clean

clean:
	rm -f  $(OBJECTS) 
	@echo "cleanup done!"

.PHONEY: remove

remove: clean
	rm -f  $(BINDIR)/main
	@echo "executable removed!"
