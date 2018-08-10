#!/bin/make

PROJ_DIR = $(HOME)/specialfunctions

PROJ_INC = $(PROJ_DIR)/inc

PROJ_SRC = $(PROJ_DIR)/src

PROJ_BIN = $(PROJ_DIR)/bin

#

GSL_DIR = $(PROJ_DIR)/third_party/gsl

GSL_INC = $(GSL_DIR)/include 

GSL_LIB = $(GSL_DIR)/lib 

#

CXX = g++

CXXOPTS = -O3 -std=c++0x -v 

#
#---- Includes 
#

INCS = -I$(PROJ_INC) -I$(GSL_INC) 

#
#---- Libraries for the linker 
#
#     NB: Some of these are shared, so will need 
#         incorporated in LD_LIBRARY_PATH
#

LIBDIRS = -L $(GSL_LIB) 

LIBS=-lm -lrt -lgsl -lgslcblas -lpthread  

#

SRCS=$(PROJ_SRC)/assoc_legendre_tests.cxx

EXEC=$(PROJ_BIN)/test.x

#
#---- Build rules 
#

all: clean test1 test2

test1:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/test1.x $(PROJ_SRC)/assoc_legendre_tests.cxx $(LIBS) 

test2:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/test2.x $(PROJ_SRC)/test2.cxx $(LIBS) 


clean:
	$(RM) -vf *.x *.o 

#
#---- End of file 
#
