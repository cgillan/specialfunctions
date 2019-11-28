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
   //
   //---- Optionally select to compute one or other kind only
   //

   bool b_compute_first_kind =   true;

   bool b_compute_second_kind = false;

   //
   //---- Optionally monitor floating point exceptions
   //

   bool b_monitor_fl_pt_exceptions = false;

   //
   //---- Banner header 
   //

   printf("\n\n");
   printf("     Computation of associated Legendre functions \n");
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
    xarg = 0.00e+00; xarg_vec.push_back(xarg);

    xarg = 0.25e+00; xarg_vec.push_back(xarg);

    xarg = 0.50e+00; xarg_vec.push_back(xarg);

    xarg = 0.75e+00; xarg_vec.push_back(xarg);

    xarg = 0.90e+00; xarg_vec.push_back(xarg);

    xarg = 1.25e+00; xarg_vec.push_back(xarg);

    xarg = 1.50e+00; xarg_vec.push_back(xarg);

    xarg = 1.75e+00; xarg_vec.push_back(xarg);

    xarg = 2.00e+00; xarg_vec.push_back(xarg);

    xarg = 2.25e+00; xarg_vec.push_back(xarg);

    xarg = 2.50e+00; xarg_vec.push_back(xarg);

    xarg = 2.75e+00; xarg_vec.push_back(xarg);

    xarg = 3.00e+00; xarg_vec.push_back(xarg);

    xarg = 3.25e+00; xarg_vec.push_back(xarg);

    xarg = 3.50e+00; xarg_vec.push_back(xarg);

    xarg = 3.75e+00; xarg_vec.push_back(xarg);

    xarg = 4.00e+00; xarg_vec.push_back(xarg);

    xarg = 4.25e+00; xarg_vec.push_back(xarg);

    xarg = 4.50e+00; xarg_vec.push_back(xarg);

    xarg = 4.75e+00; xarg_vec.push_back(xarg);

    xarg = 5.00e+00; xarg_vec.push_back(xarg);

    xarg = 5.25e+00; xarg_vec.push_back(xarg);

    xarg = 5.50e+00; xarg_vec.push_back(xarg);

    xarg = 5.75e+00; xarg_vec.push_back(xarg);

    xarg = 6.00e+00; xarg_vec.push_back(xarg);

    xarg = 6.25e+00; xarg_vec.push_back(xarg);

    xarg = 6.50e+00; xarg_vec.push_back(xarg);

    xarg = 6.75e+00; xarg_vec.push_back(xarg);

    xarg = 7.00e+00; xarg_vec.push_back(xarg);

    xarg = 7.25e+00; xarg_vec.push_back(xarg);

    xarg = 7.50e+00; xarg_vec.push_back(xarg);

    xarg = 7.75e+00; xarg_vec.push_back(xarg);

    xarg = 8.00e+00; xarg_vec.push_back(xarg);

    xarg = 8.25e+00; xarg_vec.push_back(xarg);

    xarg = 8.50e+00; xarg_vec.push_back(xarg);

    xarg = 8.75e+00; xarg_vec.push_back(xarg);

    xarg = 9.00e+00; xarg_vec.push_back(xarg);

    xarg = 9.25e+00; xarg_vec.push_back(xarg);

    xarg = 9.50e+00; xarg_vec.push_back(xarg);

    xarg = 9.75e+00; xarg_vec.push_back(xarg);

    xarg =10.00e+00; xarg_vec.push_back(xarg);
   }

   //
   //---- Set the maximum L and M   
   //

   int const Lmax = 30;
   int const Mmax = 30;

   printf("\n\n");
   printf("     Maximum L value = %4d \n", Lmax);
   printf("     Maximum M value = %4d \n", Mmax);

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
       //---- Regular functions at this argument 
       //

       if(b_compute_first_kind)
       {

       for(int m=0; m<=Mmax; ++m)
          {
           printf("     Computing regular Legendre functions P_lm(x) for m value = %d", m);
           printf("\n\n");

           //

           std::vector<long double> plm_vec;

           plm_vec.resize(Lmax+1);

           //
           //..... Compute for l=m, ...., Lmax
           //

           std::feclearexcept(FE_ALL_EXCEPT);

           unnormalised_associated_regular_Legendre(Lmax,m,x,plm_vec);

           if(b_monitor_fl_pt_exceptions)
             {
              monitor_fl_pt_exceptions();
             }

           //

           //printf("     Real argument (x) = %15.6Lf ", x);
           //printf("\n\n");
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

       }
        // End of optional selection on computing first kind of functions

       //
       //---- Irregular functions at this argument 
       //

       if(b_compute_second_kind)
         {
          std::vector<std::vector<long double> > qlm_mat;

          qlm_mat.resize(Lmax+1);

          for(int ll=0; ll<=Lmax; ++ll)
             {
              qlm_mat[ll].resize(Mmax+1);
             }
 
          //

          std::feclearexcept(FE_ALL_EXCEPT);

          unnormalised_associated_irregular_Legendre(Mmax,Lmax,x,qlm_mat);

          if(b_monitor_fl_pt_exceptions)
             {
              monitor_fl_pt_exceptions();
             }

          //

          printf("     Real argument (x) = %15.6Lf ", x);
          printf("\n\n");
          printf("     Computed associated Legendre functions of the second kind (irregular)");
          printf("\n\n");
          printf("      l    m        x         Associated Legendre function \n");
          printf("     ---  ---  -------------  ---------------------------- \n");

          for(int ll=0; ll<=Lmax; ++ll)
             {
              for(int mm=0; mm<=Mmax; ++mm)
                 {
                  printf("     %3d  %3d  %13.7Lf  %28.14Le \n", ll, mm, x, qlm_mat[ll][mm]);
                 }
             }
         }
          // End of optional selection on computing second kind of functions
      }
       // End of loop over arguments "x"

   //
   //---- End of main program 
   //

   printf("\n\n");
   printf("     **** End of computation of associated Legendre functions");
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

