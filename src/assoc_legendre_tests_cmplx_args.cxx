//************************************************************************
//************************************************************************
//
//   assoc_legendre_tests.cxx  
//
//   NB: Need --std=c++11 (or similar) on 
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
#include <string>
#include <sstream>

#include <cmath>
#include <cfenv>
#include <climits> 
#include <cfloat> 
#include <vector>
#include <algorithm>

#include "associated_legendre_functions.hxx"
  
#pragma STDC FENV_ACCESS ON  // Enable catching floating point exceptions

static void monitor_fl_pt_exceptions();

/**
 *   Main program - test harness
 *
 */

int main(int argc, char **argv)
  {
   printf("\n\n");
   printf("     Computation of associated legendre functions \n");
   printf("     --------------------------------------------   ");
   printf("\n\n");
   printf("     Real argument and integer order and degree ");
   printf("\n\n");

   //
   
   printf("     Numeric limits on this compiler and operating system");
   printf("\n\n");
   printf("     int:             \n"); 
   printf("        minimum value: %d \n", std::numeric_limits<int>::min());    
   printf("        maximum value: %d \n", std::numeric_limits<int>::max());     
   printf("        int is signed: %d \n", std::numeric_limits<int>::is_signed);  
   printf("        non-sign bits: %d \n", std::numeric_limits<int>::digits);      
   printf("        has infinity:  %d \n", std::numeric_limits<int>::has_infinity); 
   printf("     float:           \n");
   printf("         minimum value : %13.7e \n", std::numeric_limits<float>::min());
   printf("         maximum value : %13.7e \n", std::numeric_limits<float>::max());
   printf("         epsilon       : %13.7e \n", std::numeric_limits<float>::epsilon());
   printf("         exponent radix: %d     \n", FLT_RADIX);
   printf("     double:           \n");
   printf("         minimum value : %13.7e \n", std::numeric_limits<double>::min());
   printf("         maximum value : %13.7e \n", std::numeric_limits<double>::max());
   printf("         epsilon       : %13.7e \n", std::numeric_limits<double>::epsilon());
   printf("     long double:      \n");
   printf("         minimum value : %13.7Le \n", std::numeric_limits<long double>::min());
   printf("         maximum value : %13.7Le \n", std::numeric_limits<long double>::max());
   printf("         epsilon       : %13.7Le \n", std::numeric_limits<long double>::epsilon());
   printf("\n");

   //
   //---- Set the maximum L and M and argument "x"  
   //

   int const Lmax =  200;
   int const Mmax =  200;

   long double const x = 1.10e+00;

   //======================================================================
   //
   //     R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================

   for(int m=0; m<=Mmax; ++m)
      {
       printf("\n\n");
       printf("     Computing regular associated Legendre functions \n");
       printf("     for m = %d and l in the range [%d,%d] ", m,m,Lmax);

       //

       std::vector<long double> plm_vec;

       plm_vec.resize(Lmax+1);

       //
       //..... Compute for l=m, ...., Lmax
       //

       std::feclearexcept(FE_ALL_EXCEPT);

       unnormalised_associated_regular_Legendre(Lmax,m,x,plm_vec);

       monitor_fl_pt_exceptions();

       //

       printf("\n\n");
       printf("     Real argument (x) = %15.6Lf ", x);
       printf("\n\n");
       printf("     Computed associated Legendre functions of the first kind (regular)");
       printf("\n\n");
       printf("      l    m        x        Associated Legendre function \n");
       printf("     ---  ---  ------------- ---------------------------- \n");

       for(int l=m; l<=Lmax; ++l)
          {
           printf("     %3d  %3d  %15.6Le       %15.6Le \n", l, m, x, plm_vec[l]);
          }
      }
       // End loop over m values 

   //======================================================================
   //
   //     I R R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================

   printf("\n\n");
   printf("     Computing regular associated Legendre functions \n");
   printf("     for m in range [%d,%d] and l in the range [%d,%d] ", 
           0,Mmax, 0,Lmax);

   //
   //---- Create rectangular storage
   //

   std::vector<std::vector<long double> > qlm_mat;

   qlm_mat.resize(Lmax+1);

   for(int ll=0; ll<=Lmax; ++ll)
      {
       qlm_mat[ll].resize(Mmax+1);
      }
 
   //

   std::feclearexcept(FE_ALL_EXCEPT);

   unnormalised_associated_irregular_Legendre_big_arg(Lmax,Mmax,x,qlm_mat);

   monitor_fl_pt_exceptions();

   //

   printf("\n\n");
   printf("     Real argument (x) = %15.6Lf ", x);
   printf("\n\n");
   printf("     Computed associated Legendre functions of the second kind (irregular)");
   printf("\n\n");
   printf("      l    m        x        Associated Legendre function \n");
   printf("     ---  ---  ------------- ---------------------------- \n");

   for(int ll=0; ll<=Lmax; ++ll)
      {
       for(int mm=0; mm<=Mmax; ++mm)
          {
           printf("     %3d  %3d  %15.6Le       %15.6Le \n", ll, mm, x, qlm_mat[ll][mm]);
          }
      }

   //
   //---- End of main program 
   //

   printf("\n\n");
   printf("     **** End of computation of associated legendre functions");
   printf("\n\n");

   return 0;
  }
   // End of main

/**
 *   monitor_fl_pt_exceptions()
 *
 *   Test current state of the floating point flags and report
 *   to the job log
 * 
 *   Before calling the code to be monitored, all floating point 
 *   exceptions whould be cleared by calling 
 *
 *     std::feclearexcept(FE_ALL_EXCEPT);
 *
 *   This routine can be called after the critial code to 
 *   report on what happened.  
 *
 */
static void monitor_fl_pt_exceptions()
  {
   printf("\n\n");
   printf("     Reporting status of floating point exception flags.");
   printf("\n\n");

   if(std::fetestexcept(FE_ALL_EXCEPT))
     {
      if(std::fetestexcept(FE_DIVBYZERO)) 
        {
         printf("     **** Division by zero reported");
         printf("\n");
        }

      if(std::fetestexcept(FE_INEXACT)) 
        {
         printf("     **** Inexact computation ");
         printf("\n");
        }

      if(std::fetestexcept(FE_INVALID)) 
        {
         printf("     **** Invalid computation ");
         printf("\n");
        }

      if(std::fetestexcept(FE_OVERFLOW)) 
        {
         printf("     **** Overflow reported ");
         printf("\n");
        }

      if(std::fetestexcept(FE_UNDERFLOW)) 
        {
         printf("     **** Underflow reported ");
         printf("\n");
        }
     }
   else
     {
      printf("     No floating point exception flags are set.");
      printf("\n");
     }
      // End of test on floating point exceptions being raised

   //

   return;
  }
   // End of monitor_fl_pt_exceptions()

//************************************************************************
//************************************************************************
//
//   End of file 
//
//************************************************************************
//************************************************************************

