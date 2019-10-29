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
#include <typeinfo>

#include <cmath>
#include <cfenv>
#include <climits> 
#include <cfloat> 
#include <vector>
#include <algorithm>

#include <associated_legendre_functions_cmplx.hxx>  

#include <monitor_fl_pt_exceptions.hxx>

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
   printf("     Complex argument and integer order and degree ");
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

   std::complex<long double> zarg;
 
   zarg.real(2.00e+00); zarg.imag(0.0e+00);

   //======================================================================
   //
   //     R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================

   std::vector<std::vector<std::complex<long double> > > plm_mat;

   plm_mat.resize(Lmax+1);

   for(int l=0; l<=Lmax; ++l)
      {
       plm_mat[l].resize(l+1);
      }

   std::feclearexcept(FE_ALL_EXCEPT);

   complex_unnormalized_assoc_regular_legendre((unsigned int) Lmax,
                                               (unsigned int) Mmax,
                                               zarg,
                                               plm_mat);

   monitor_fl_pt_exceptions();

   //

   printf("\n\n");
   printf("     Computed associated Legendre functions of the first kind (regular)");
   printf("\n\n");
   printf("      l    m           Argument (z)                  Associated Legendre function \n");
   printf("     ---  ---  ---------------------------------  --------------------------------- \n");

   std::string cformat_str = " ";

   for(int l=0; l<=Lmax; ++l)
      {
       for(int m=0; m<=l; ++m)
          {
           long double const xarg = zarg.real();
           long double const yarg = zarg.imag();

           //

           std::complex<long double> const zplm = plm_mat[l][m];

           long double const xplm = zplm.real();
           long double const yplm = zplm.imag();

           printf("     %3d  %3d  (%15.6Le,%15.6Le)   (%15.6Le, %15.6Le) \n", 
                   l, m, 
                   xarg, yarg, 
                   xplm, yplm);
          }
      }

   //======================================================================
   //
   //     I R R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================

   printf("\n\n");
   printf("     Computing regular associated Legendre functions \n");
   printf("     for m in range [%d,%d] and l in the range [%d,%d] ", 0,Mmax, 0,Lmax);

   //

   std::vector<std::vector<std::complex<long double> > > qlm_mat;

   qlm_mat.resize(Lmax+1);

   for(int ll=0; ll<=Lmax; ++ll)
      {
       qlm_mat[ll].resize(Mmax+1);
      }
 
   //

   std::feclearexcept(FE_ALL_EXCEPT);

   complex_unnormalized_assoc_irregular_legendre(Lmax,Mmax,zarg,qlm_mat);

   monitor_fl_pt_exceptions();

   //

   printf("\n\n");
   printf("     Complex argument (z) = (%15.6Lf,%15.6Lf) ", zarg.real(), zarg.imag());
   printf("\n\n");
   printf("     Computed associated Legendre functions of the second kind (irregular)");
   printf("\n\n");
   printf("      l    m        z        Associated Legendre function \n");
   printf("     ---  ---  ------------- ---------------------------- \n");

   for(int ll=0; ll<=Lmax; ++ll)
      {
       for(int mm=0; mm<=Mmax; ++mm)
          {
           printf("     %3d  %3d  (%15.6Le,%15.6Le)   (%15.6Le,%15.6Le) \n", 
                   ll, mm, 
                   zarg.real(), 
                   zarg.imag(),
                   qlm_mat[ll][mm].real(), 
                   qlm_mat[ll][mm].imag());
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

//************************************************************************
//************************************************************************
//
//   End of file 
//
//************************************************************************
//************************************************************************

