TOOLBOXDIR = ./toolbox

CCFILES =          compute.cpp

CFILES =           main.c \
				   create_variables.c \
				   permeabilities.c \
				   set_boundaries.c \
				   set_source_term.c \
				   solve_linearsystem_AMG.c \
				   flux_reconstruction.c \
				   free_variables.c \
				   posprocessing_output.c \
                                   karst.c\
				   vtk.c \
				   reservoir.c

EXE = exe 
EDIR = ./run
ODIR = ./obj


CXX        := g++
CC         := gcc
LINK       := g++ -fPIC

CXXFLAGS = -O1 -funroll-all-loops -DNOZEROCOPY
CFLAGS = -O1 -g -W
AR = ar cr

# Compile files

C_DEPS = structs.h \
	     function.h \
		 variables.h

OBJS +=  $(patsubst %.cpp,$(ODIR)/%.cpp.o,$(notdir $(CCFILES)))
OBJS +=  $(patsubst %.c,$(ODIR)/%.c.o,$(notdir $(CFILES)))

$(ODIR)/%.c.o : %.c 
	$(CXX) $(CFLAGS) -o $@ -c $<

$(ODIR)/%.cpp.o : %.cpp 
	$(CXX) $(CXXFLAGS) -o $@ -c $<

default: $(OBJS)
	$(LINK) -o $(EDIR)/$(EXE) $(OBJS)


clean:
	rm -f $(ODIR)/*.c.o
	rm -f $(ODIR)/*.cpp.o

