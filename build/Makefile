shell := /bin/bash

# Use this utility to perform a make from the command-line using mpicc
GRIDBASE  ?= $(shell grid-config --prefix)
GSLBASE=/home/dp008/dp008/dc-kett1/install/gsl

# grid-confic --cxx not working so hard code same compilers as grid
CC=icc 
CXX=icpc 
MPICXX=mpiicpci


#CXX       := $(shell grid-config --cxx)
CXXFLAGS  := $(shell grid-config --cxxflags) -I$(GRIDBASE)/include -I$(GSLBASE)/include 
LDLIBS    := $(shell grid-config --libs)    -lGrid -lgsl 
LDFLAGS   := $(shell grid-config --ldflags) -L$(GRIDBASE)/lib -L$(GSLBASE)/lib

# set where source files are
SRC_DIRS ?= ../src
VPATH :=$(SRC_DIRS)

# define sources, objects and progs
SRCFILES := $(shell find $(SRC_DIRS) -name *.cc )
OBJFILES := $(notdir $(patsubst %.cc, %.o, $(SRCFILES)))
PROGFILES:= $(notdir $(patsubst %.cc, %, $(SRCFILES)))


INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

all: $(PROGFILES)

.PHONY: all clean

all: $(PROGFILES)
clean:
	 $(RM) $(OBJFILES) $(PROGFILES) *~









