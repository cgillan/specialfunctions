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

PROJ_LIB = $(PROJ_DIR)/lib

PROJ_BIN = $(PROJ_DIR)/bin

#

CXX = g++

CXXOPTS = -O3 -std=c++0x 

CXXOPTS += -v 

#

AR = ar 

AROPTS = -crv 

#
#---- Includes 
#

INCS = -I$(PROJ_INC)  

#
#---- Libraries for the linker 
#

LIBDIRS = -L $(PROJ_LIB)

LIBS= -lmonitor -lm -lrt 

#
#---- Build rules 
#
#     Separate rule for each teet harness
#

all: clean monitorlib cplm testreal testcmplx 

monitorlib:
	$(CXX) $(INCS) $(CXXOPTS) -c   \
               -o $(PROJ_SRC)/monitor_fl_pt_exceptions.o \
               $(PROJ_SRC)/monitor_fl_pt_exceptions.cxx 
	$(AR)  $(AROPTS) $(PROJ_LIB)/libmonitor.a \
               $(PROJ_SRC)/monitor_fl_pt_exceptions.o    

cplm:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS)  \
               -o $(PROJ_BIN)/cplm.x \
               $(PROJ_SRC)/cplm.cxx $(LIBS) 

testreal:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/assoc_legendre_tests_real_args.x \
                  $(PROJ_SRC)/assoc_legendre_tests_real_args.cxx $(LIBS) 

testcmplx:
	$(CXX) $(INCS) $(CXXOPTS) $(LIBDIRS) \
               -o $(PROJ_BIN)/assoc_legendre_tests_cmplx_args.x \
                  $(PROJ_SRC)/assoc_legendre_tests_cmplx_args.cxx $(LIBS) 

clean:
	$(RM) -vf $(PROJ_BIN)/*.x $(PROJ_SRC)/*.o $(PROJ_LIB)/*.a 

#
#---- End of file 
#
