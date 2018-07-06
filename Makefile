#!/bin/make

PROJ_DIR = $(HOME)/CPC_paper_proposal

PROJ_INC = $(PROJ_DIR)/inc

PROJ_SRC = $(PROJ_DIR)/src

PROJ_BIN = $(PROJ_DIR)/bin

#

CXX = g++

CXXOPTS = -O2   

CXXOPTS +=-std=c++0x -v 

LIBS = -lm 
#LIBS += -lrt   

#
#---- Files to compile and link to binaries
#

SRCS=$(PROJ_SRC)/assoc_legendre_tests.cxx

EXEC=$(PROJ_BIN)/test.x

#
#---- Build rules 
#

all: clean complink 

complink:
	$(CXX) -I $(PROJ_INC) $(CXXOPTS) -o $(EXEC) $(SRCS) $(LIBS) 

clean:
	$(RM) -vf *.x *.o 

#
#---- End of file 
#
