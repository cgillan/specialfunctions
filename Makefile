#!/bin/make
#
#  Note the assumption here is that the present working 
#  directory is at the top of the project sub-directory 
#  hierarchy.
#
#  Source code is located in ./src
#
#  Exdcutables are placed in ./bin
#
#

PROJ_DIR = $(PWD)

PROJ_INC = $(PROJ_DIR)/inc

PROJ_SRC = $(PROJ_DIR)/src

PROJ_BIN = $(PROJ_DIR)/bin

#

CXX = g++

CXXOPTS = -O3 -std=c++0x 

CXXOPTS += -v 

#
#---- Includes 
#

INCS = -I$(PROJ_INC)  

#
#---- Libraries for the linker 
#

LIBS=-lm -lrt 

#
#---- Build rules 
#
#     Separate rule for each teet harness
#

all: clean testzhangjin testreal testcmplx 

testzhangjin:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS)  \
               -o $(PROJ_BIN)/zhang_jin_plm.x \
               $(PROJ_SRC)/zhang_jin_plm.cxx $(LIBS) 

testreal:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/assoc_legendre_tests_real_args.x \
                  $(PROJ_SRC)/assoc_legendre_tests_real_args.cxx $(LIBS) 

testcmplx:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/assoc_legendre_tests_cmplx_args.x \
                  $(PROJ_SRC)/assoc_legendre_tests_cmplx_args.cxx $(LIBS) 

clean:
	$(RM) -vf $(PROJ_BIN)/*.x $(PROJ_SRC)/*.o 

#
#---- End of file 
#
