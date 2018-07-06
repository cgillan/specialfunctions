#!/bin/make

PROJ_DIR = .

PROJ_INC = $(PROJ_DIR)/inc

PROJ_SRC = $(PROJ_DIR)/src

PROJ_BIN = $(PROJ_DIR)/bin

#

CXX = g++

CXXOPTS = -O2

CXXOPTS +=-std=c++0x -g3

LIBS = -lm
#LIBS += -lrt


#
#---- Build rules
#

all: bin/test.x bin/benchmark.x

bin/test.x: $(PROJ_SRC)/assoc_legendre_tests.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) $(PROJ_SRC)/assoc_legendre_tests.cxx -o bin/test.x

bin/benchmark.x: $(PROJ_SRC)/assoc_legendre_benchmark.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) $(PROJ_SRC)/assoc_legendre_benchmark.cxx -o bin/benchmark.x


.PHONY: clean
clean:
	$(RM) -vrf bin
	$(RM) -vrf prof

.PHONY: run
run: bin/test.x
	./bin/test.x



#
#---- Profiling rules
#

.PHONY: prof
prof: bin/benchmark.x
	mkdir -p prof
	LD_PRELOAD=/usr/lib/libprofiler.so CPUPROFILE=prof/cpu_profile bin/benchmark.x
	google-pprof --text bin/benchmark.x prof/cpu_profile


#
#---- End of file
#
