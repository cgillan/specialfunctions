//************************************************************************
//************************************************************************
//
//   assoc_legendre_tests.cxx  
//
//   NB: Need --std=c++0x (or similar) on 
//       g++ otherwise this fails 
//
//   Copyright (c) 2018 Charles J Gillan  
//   All rights reserved
//
//************************************************************************
//************************************************************************

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits> 

#include <string>

#include <vector>
#include <algorithm>

#include "associated_legendre_functions.hxx"
  
/**
 *   Main program - test harness
 *
 */

int main(int argc, char **argv)
  {
   int const Lmax =  6;
   int const Mmax =  6;

   //
   //---- Set the argument 
   //

   double const x = 5.0e+00;

   //======================================================================
   //
   //     R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================
   //
   //---- Create lower triangular storage
   //

   std::vector<std::vector<double> > plm_lower_triangle_mat;

   plm_lower_triangle_mat.resize(Lmax+1);

   for(int ll=0; ll<=Lmax; ++ll)
      {
       plm_lower_triangle_mat[ll].resize(ll+1);
      }
 
   //
   //---- Use the recurrence relationship to build up P_lm 
   //

   for(int m=0; m<=Mmax; ++m)
      {
       std::vector<double> plm_vec;

       plm_vec.resize(Lmax+1);

       //
       //..... Compute for l=m, ...., Lmax
       //

       unnormalised_associated_regular_Legendre(Lmax,m,x,plm_vec);

       //
       //..... Copy to triangle storage
       //

       for(int l=m; l<=Lmax; ++l)
          {
           plm_lower_triangle_mat[l][m] = plm_vec[l];
          }
      }

   //
   //---- Print lower triangular matrix
   //

   std::string const header_str = "P_lm values";

   std::cout << "\n\n";
   std::cout << " -------------------------------------------------- ";
   std::cout << "\n\n";
   std::cout << "     Matrix of P_{lm} values for argument (x) = " << x; 
   std::cout << "\n\n";
   std::cout << "     Each column is an m value, each row an l value";

   printLowerTriangularPlmMat(plm_lower_triangle_mat, Lmax);
   
   //======================================================================
   //
   //     I R R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================
   //
   //---- Create rectangular storage
   //

   std::vector<std::vector<double> > qlm_mat;

   qlm_mat.resize(Lmax+1);

   for(int ll=0; ll<=Lmax; ++ll)
      {
       qlm_mat[ll].resize(Mmax+1);
      }
 
   //

   unnormalised_associated_irregular_Legendre_big_arg(Lmax,Mmax,x,qlm_mat);

   //

   std::cout << "\n\n";
   std::cout << " -------------------------------------------------- ";
   std::cout << "\n\n";
   std::cout << "     Matrix of Q_{lm} values for argument (x) = " << x; 
   std::cout << "\n\n";

   printRectangularQlmMat(Lmax,Mmax,qlm_mat);
  }
   // End of main

//************************************************************************
//************************************************************************
//
//   End of file 
//
//************************************************************************
//************************************************************************

