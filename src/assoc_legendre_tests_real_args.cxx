//************************************************************************
//************************************************************************
//
//   assoc_legendre_tests_real_args.cxx  
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

#include "monitor_fl_pt_exceptions.hxx"

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
   //---- Determine, at runtime, the data type for which we were compiled
   //

   long double xarg;

   std::string cformat_type_str = " ";

   if(typeid(xarg) == typeid( float ))
     {
      cformat_type_str = "     Data type for argument is: float. \n\n     File name = %s";
     }
   else if(typeid(xarg) == typeid( double ))
     {
      cformat_type_str = "     Data type for argument is: double. \n\n     File name = %s";
     }
   else if(typeid(xarg) == typeid( long double ))
     {
      cformat_type_str = "     Data type for argument is: long double. \n\n     File name = %s";
     }
   else
     {
      printf("\n\n");
      printf("     **** Error: Type for the complex argument is unknown");
      printf("\n\n");

      exit(0);
     }

   printf(cformat_type_str.c_str(), __FILE__);
   printf("\n\n");

   //
   //---- Prepare vector of arguments 
   //

   std::vector<long double> xarg_vec;

   {
    xarg = 0.5e+00; xarg_vec.push_back(xarg);

    xarg = 1.5e+00; xarg_vec.push_back(xarg);

    xarg = 2.0e+00; xarg_vec.push_back(xarg);
   }

   //
   //---- Set the maximum L and M   
   //

   int const Lmax = 20;
   int const Mmax = 20;

   printf("\n\n");
   printf("      Maximum L value = %4d \n", Lmax);
   printf("      Maximum M value = %4d \n", Mmax);

   //================================================================
   //
   //    L O O P  O V E R   A R G U M E N T S
   //
   //================================================================

   for(int indx=0; indx<xarg_vec.size(); ++indx)
      {
       auto const x = xarg_vec[indx];

       printf("\n\n "); 
          for(int icol=2; icol<72;++icol) printf("-"); 

       printf("\n\n");

       //
       //---- Regular functions 
       //

       for(int m=0; m<=Mmax; ++m)
          {
           printf("     Computing regular Legendre functions P_lm(x) for m value = %d", m);
           printf("\n\n");
           printf("     Argument (x) = %15.6Lf", x);
           printf("\n\n");

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
           printf("      l    m        x         Associated Legendre function \n");
           printf("     ---  ---  -------------  ---------------------------- \n");

           for(int l=m; l<=Lmax; ++l)
              {
               printf("     %3d  %3d  %13.7Lf  %28.14Le \n", l, m, x, plm_vec[l]);
              }

           printf("\n\n");
          }
           // End loop over m values 

       //
       //---- Irregular functions 
       //

       std::vector<std::vector<long double> > qlm_mat;

       qlm_mat.resize(Lmax+1);

       for(int ll=0; ll<=Lmax; ++ll)
          {
           qlm_mat[ll].resize(Mmax+1);
          }
 
       //

       std::feclearexcept(FE_ALL_EXCEPT);

       //unnormalised_associated_irregular_Legendre_big_arg(Mmax,Lmax,x,qlm_mat);

       real_unnormalized_assoc_irregular_legendre(Mmax,Lmax,x,qlm_mat);

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
      }
       // End of loop over arguments "x"

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

