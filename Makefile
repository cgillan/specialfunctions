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

all: bin/benchmark_seq.x bin/benchmark_omp.x bin/benchmark_with_results.x

bin/benchmark_seq.x: Makefile $(PROJ_SRC)/benchmark.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) $(PROJ_SRC)/benchmark.cxx -o bin/benchmark_seq.x -DSKIP_PRINTING_RESULTS=1

bin/benchmark_omp.x: Makefile $(PROJ_SRC)/benchmark.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) -fopenmp $(PROJ_SRC)/benchmark.cxx -o bin/benchmark_omp.x -DSKIP_PRINTING_RESULTS=1

bin/benchmark_with_results.x: Makefile $(PROJ_SRC)/benchmark.cxx $(PROJ_INC)/associated_legendre_functions.hxx
	mkdir -p bin
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) -fopenmp $(PROJ_SRC)/benchmark.cxx -o bin/benchmark_with_results.x


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

.PHONY: prof_seq
prof_seq: bin/benchmark_seq.x
	mkdir -p prof
	LD_PRELOAD=/usr/lib/libprofiler.so CPUPROFILE=prof/cpu_profile bin/benchmark_seq.x
	google-pprof --text bin/benchmark_seq.x prof/cpu_profile
	google-pprof --web bin/benchmark_seq.x prof/cpu_profile

.PHONY: prof_omp
prof_omp: bin/benchmark_omp.x
	mkdir -p prof
	LD_PRELOAD=/usr/lib/libprofiler.so CPUPROFILE=prof/cpu_profile bin/benchmark_omp.x
	google-pprof --text bin/benchmark_omp.x prof/cpu_profile
	google-pprof --web bin/benchmark_omp.x prof/cpu_profile


#
#---- End of file
#
