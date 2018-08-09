#!/bin/make

PROJ_DIR = .

PROJ_INC = $(PROJ_DIR)/inc

PROJ_SRC = $(PROJ_DIR)/src

PROJ_BIN = $(PROJ_DIR)/bin

#

CXX = g++

CXXOPTS = -O3 -std=c++14 -g3

LIBS = -lm
#LIBS += -lrt


#
#---- Build rules
#

all: bin/benchmark.x bin/benchmark_with_results.x

bin/benchmark.x: Makefile $(PROJ_SRC)/benchmark.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) $(PROJ_SRC)/benchmark.cxx -o bin/benchmark.x -DSKIP_PRINTING_RESULTS=1

bin/benchmark_with_results.x: Makefile $(PROJ_SRC)/benchmark.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) $(PROJ_SRC)/benchmark.cxx -o bin/benchmark_with_results.x


.PHONY: clean
clean:
	$(RM) -vrf bin
	$(RM) -vrf prof

.PHONY: run
run: bin/benchmark_with_results.x
	./bin/benchmark_with_results.x



#
#---- Profiling rules
#

.PHONY: prof
prof: bin/benchmark.x
	mkdir -p prof
	LD_PRELOAD=/usr/lib/libprofiler.so CPUPROFILE=prof/cpu_profile bin/benchmark.x
	google-pprof --text bin/benchmark.x prof/cpu_profile
	google-pprof --web bin/benchmark.x prof/cpu_profile


#
#---- End of file
#
