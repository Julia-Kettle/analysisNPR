# Use this utility to perform a make from the command-line using mpicc
GRIDBASE  ?= $(shell grid-config --prefix)

# grid-confic --cxx not working so hard code same compilers as grid
CC=icc 
CXX=icpc 
MPICXX=mpiicpci


#CXX       := $(shell grid-config --cxx)
CXXFLAGS  := $(shell grid-config --cxxflags) -I$(GRIDBASE)/include 
LDLIBS    := $(shell grid-config --libs)    -lGrid 
LDFLAGS   := $(shell grid-config --ldflags) -L$(GRIDBASE)/lib 

TARGET ?= bilinearAnalysis
SRC_DIRS ?= ../src

SRCS := $(shell find $(SRC_DIRS) -name *.cc -or -name *.c -or -name *.s)
OBJS := $(addsuffix .o,$(basename $(SRCS)))
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

$(TARGET): $(OBJS)
		$(CC) $(LDFLAGS) $(OBJS) -o $@ $(LOADLIBES) $(LDLIBS)

.PHONY: clean
	clean:
		$(RM) $(TARGET) $(OBJS) $(DEPS)

-include $(DEPS)








