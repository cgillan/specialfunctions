#!/bin/make

PROJ_DIR = $(PWD)

PROJ_INC = $(PROJ_DIR)/inc

PROJ_SRC = $(PROJ_DIR)/src

PROJ_BIN = $(PROJ_DIR)/bin

#

CXX = g++

CXXOPTS = -O3 -std=c++0x  

#
#---- Includes 
#

INCS = -I$(PROJ_INC) 

#
#---- Libraries for the linker 
#
#     NB: Some of these are shared, so will need 
#         incorporated in LD_LIBRARY_PATH
#

LIBDIRS = -L . 

LIBS=-lm -lrt -lpthread  

#
#---- Build rules 
#

all: clean testreal testcmplx

testreal:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/testreal.x     \
                  $(PROJ_SRC)/assoc_legendre_tests_real_args.cxx $(LIBS) 

testcmplx:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/testcmplx.x    \
                  $(PROJ_SRC)/assoc_legendre_tests_cmplx_args.cxx $(LIBS) 


clean:
	$(RM) -vf $(PROJ_BIN)/*.x *.o 

#
#---- End of file 
#
