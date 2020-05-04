# specialfunctions

Background
----------

This directory hierarchy contains the code which implements 
computation of the associated Legendre functions for integer
only values of the degree and order, which are typically denoted
as "l" and "m" respectively. 

The code is released under the GPL3 license. 

Directory structure
-------------------

./inc  The C++ templates as described above are located here

./src  A test driver programs which use the templates are placed here

       This driver computes the P_lm and Q_lm and may therefore
       be used to compare with the table of function values 
       reported in the research paper above.

./bin  The binary for the test programs are placed in this location

Compilation
-----------

There is a Makefile at the top level of this directory structure.

All that is required is to change to that directory and to type 

      make
      
Note that the C++11 features in the compiler must be enabled.

The code has been tested with g++ and so the --std=c++0x option 
needs to be used on the compiler line.
The test binary can then be run by typing, for example 

      ./bin/test2.x 

References to related work
--------------------------

 1.   A new Fortran 90 program to compute regular and irregular
      associated Legendre functions
      B L Schneider, J Segura, A Gil, X Guan and K Bartschat
      Journal Computer Physics  Comminucations
      Volume 181, 2010, pp 2091-7

 2.   Computation of Special Functions 
      Shanjie Zhang aand Jianming Jin
      pub: John Wiley & Sons Inc., 1996,
      ISBN 0-471-11963-6
      QA351.C45

 3.   Study of Two Center Integrals Useful in Calculations on 
      Molecular Structure
      V. General methods for Diatomic Integrals Applicable 
         to Digital Computers  
      Journal of Chemical Physics, Vol 41., Number 9, 1964, pp 2578-99


*** End ***

