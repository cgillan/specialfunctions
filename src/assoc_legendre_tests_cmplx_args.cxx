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

   std::complex<long double> zarg;
 
   zarg.real(2.00e+00); zarg.imag(0.0e+00);

   //

   std::string cformat_type_str = " ";

   if(typeid(zarg) == typeid( std::complex<float>(0.0,1.0) ) )
     {
      cformat_type_str = "     Data type for complex argument is: float. \n\n     File name = %s";
     }
   else if(typeid(zarg) == typeid( std::complex<double>(0.0,1.0) ) )
     {
      cformat_type_str = "     Data type for complex argument is: double. \n\n     File name = %s";
     }
   else if(typeid(zarg) == typeid( std::complex<long double>(0.0,1.0) ) )
     {
      cformat_type_str = "     Data type for complex argument is: long double. \n\n     File name = %s";
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

   //======================================================================
   //
   //     R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //     T E S T E D   V A L U E S 
   //
   //======================================================================
   //
   //---- Begin scoping for P(l,m) tests
   //

   {

   //---- Following are selected values from the tables of Zhang and Jin
   //
   //     See pages 118 and following of their book.
   //

   std::vector<int> l_vec_test_p { 1, 2, 3, 10, 2, 3, 4, 10, 3, 4, 5, 10, 4, 5, 6 };
   std::vector<int> m_vec_test_p { 1, 1, 1,  1, 2, 3, 2,  2, 3, 3, 3,  3, 4, 4, 4 };

   //
   //---- Set the maximum L and M and argument "z"  
   //

   auto iter_max_l = std::max_element(std::begin(l_vec_test_p), std::end(l_vec_test_p)); 
   auto iter_max_m = std::max_element(std::begin(m_vec_test_p), std::end(m_vec_test_p)); 

   int const Ltemp = *iter_max_l;
   int const Mtemp = *iter_max_m;

   printf("     List of P(l,m) values to be tested at each argument");
   printf("\n\n");
   printf("     Index    l      m   \n");
   printf("     -----  -----  ----- \n");

   for(int i=0; i<l_vec_test_p.size(); ++i)
      {
       printf("     %5d  %5d  %5d  \n", i+1, l_vec_test_p[i], m_vec_test_p[i]);
      }
   
   printf("\n\n");
   printf("      Maximum L value in test list = %4d \n", Ltemp);
   printf("      Maximum M value              = %4d \n", Mtemp);
   printf("\n\n");

   //
   //---- Given the way that the routines work with triangles and 
   //     matrices, we find the largest of the two l,m and dimension
   //     with that 
   //

   int const Itemp = std::max(Ltemp,Mtemp);

   int const Lmax = Itemp;
   int const Mmax = Itemp;

   //

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
   //---- Print out the results 
   //

   std::string cformat_plm_str = " ";

   if(typeid(zarg) == typeid( std::complex<float>(0.0,1.0) ) )
     {
      cformat_plm_str = "     %3d  %3d  (%15.6e,%15.6e)   (%15.6e, %15.6e) \n"; 
     }
   else if(typeid(zarg) == typeid( std::complex<double>(0.0,1.0) ) )
     {
      cformat_plm_str = "     %3d  %3d  (%15.6e,%15.6e)   (%15.6e, %15.6e) \n"; 
     }
   else if(typeid(zarg) == typeid( std::complex<long double>(0.0,1.0) ) )
     {
      cformat_plm_str = "     %3d  %3d  (%15.6Le,%15.6Le)   (%15.6Le, %15.6Le) \n"; 
     }
   else
     {
      printf("\n\n");
      printf("     **** Error: Type for the complex argument is unknown");
      printf("\n\n");
     }

   printf("\n\n");
   printf("     Computed associated Legendre functions of the first kind (regular)");
   printf("\n\n");
   printf("      l    m           Argument (z)                   Associated Legendre function    \n");
   printf("     ---  ---  ---------------------------------   ---------------------------------- \n");

   for(int i=0; i<l_vec_test_p.size(); ++i)
      {
       int const l = l_vec_test_p[i];
       int const m = m_vec_test_p[i];

       if(m > l) continue; // Q_lm can have m > l but not P_lm.

       long double const xarg = zarg.real();
       long double const yarg = zarg.imag();

       //

       std::complex<long double> const zplm = plm_mat[l][m];

       long double const xplm = zplm.real();
       long double const yplm = zplm.imag();

       printf(cformat_plm_str.c_str(), l, m, xarg, yarg, xplm, yplm);
      }

   //

   }  

   // End of  scoping for P(l,m) tests

   //======================================================================
   //
   //     I R R E G U L A R   L E G E N D R E   F U N C T I O N S
   //
   //======================================================================
   //
   //---- Begin scoping for Q(l,m) tests
   //

   {

   //---- Following are selected values from the tables of Zhang and Jin
   //
   //     See pages 118 and following of their book.
   //

   std::vector<int> l_vec_test_q {  0, 1, 2, 10, 0, 1, 2, 10, 0, 1, 2, 10, 0, 1, 2, 10 };
   std::vector<int> m_vec_test_q {  1, 1, 1,  1, 2, 2, 2,  2, 3, 3, 3,  3, 4, 4, 4,  4 };

   //
   //---- Set the maximum L and M and argument "z"  
   //

   auto iter_max_l = std::max_element(std::begin(l_vec_test_q), std::end(l_vec_test_q)); 
   auto iter_max_m = std::max_element(std::begin(m_vec_test_q), std::end(m_vec_test_q)); 

   int const Ltemp = *iter_max_l;
   int const Mtemp = *iter_max_m;

   printf("     List of Q(l,m) values to be tested at each argument");
   printf("\n\n");
   printf("     Index    l      m   \n");
   printf("     -----  -----  ----- \n");

   for(int i=0; i<l_vec_test_q.size(); ++i)
      {
       printf("     %5d  %5d  %5d  \n", i+1, l_vec_test_q[i], m_vec_test_q[i]);
      }
   
   printf("\n\n");
   printf("      Maximum L value in test list = %4d \n", Ltemp);
   printf("      Maximum M value              = %4d \n", Mtemp);
   printf("\n\n");

   //
   //---- Given the way that the routines work with triangles and 
   //     matrices, we find the largest of the two l,m and dimension
   //     with that 
   //

   int const Itemp = std::max(Ltemp,Mtemp);

   int const Lmax = Itemp;
   int const Mmax = Itemp;

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
   //---- Print out the results 
   //

   std::string cformat_qlm_str = " ";

   if(typeid(zarg) == typeid( std::complex<float>(0.0,1.0) ) )
     {
      cformat_qlm_str = "     %3d  %3d  (%15.6e,%15.6e)   (%15.6e, %15.6e) \n"; 
     }
   else if(typeid(zarg) == typeid( std::complex<double>(0.0,1.0) ) )
     {
      cformat_qlm_str = "     %3d  %3d  (%15.6e,%15.6e)   (%15.6e, %15.6e) \n"; 
     }
   else if(typeid(zarg) == typeid( std::complex<long double>(0.0,1.0) ) )
     {
      cformat_qlm_str = "     %3d  %3d  (%15.6Le,%15.6Le)   (%15.6Le, %15.6Le) \n"; 
     }
   else
     {
      printf("\n\n");
      printf("     **** Error: Type for the complex argument is unknown");
      printf("\n\n");
     }

   printf("\n\n");
   printf("     Computed associated Legendre functions of the second kind (irregular)");
   printf("\n\n");
   printf("      l    m           Argument (z)                   Associated Legendre function    \n");
   printf("     ---  ---  ---------------------------------   ---------------------------------- \n");

   for(int i=0; i<l_vec_test_q.size(); ++i)
      {
       int const ll = l_vec_test_q[i];
       int const mm = m_vec_test_q[i];

       long double const xarg = zarg.real();
       long double const yarg = zarg.imag();

       //

       std::complex<long double> const zqlm = qlm_mat[ll][mm];

       long double const xqlm = zqlm.real();
       long double const yqlm = zqlm.imag();

       //

       printf(cformat_qlm_str.c_str(), ll, mm, xarg, yarg, xqlm, yqlm);
      }

   //

   }  

   // End of  scoping for Q(l,m) tests
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

