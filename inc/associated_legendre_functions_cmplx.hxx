//************************************************************************
//************************************************************************
//
//   File: associated_legendre_functions_complex.hxx  
//
//   This file implements the equations for the computation of
//   associated Legendre polynamicals of the first and second 
//   kind of integer order (l) and degree (m) and complex 
//   argument (z), as defined in the book 
//
//        Shanjie Zhang, Jianming Jin,
//        Computation of Special Functions,
//        Wiley, 1996,
//        ISBN: 0-471-11963-6,
//        LC: QA351.C45.
//
//   The output has been verfied against the tables printed for 
//   selected real values of the argument in the book.
//
//   NB: Need -std=c++11 (or similar) on g++ otherwise this fails 
//
//   Copyright (c) 2019 Charles J Gillan  
//   All rights reserved
//
//************************************************************************
//************************************************************************

#ifndef _ASSOC_LEGENDRE_TEMPLATED_COMPLEX_INCLUDE_176496_H_

#define _ASSOC_LEGENDRE_TEMPLATED_COMPLEX_INCLUDE_176496_H_  1

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <cmath>
#include <climits> 
#include <cfloat> 
#include <vector>
#include <algorithm>
#include <complex>
#include <cstdint>
#include <typeinfo>

/**
 *   real_factorial()
 *
 *   Compute factorial of "n" using floating point 
 *
 */ 
template <typename T> 
    typename std::enable_if<std::is_floating_point<T>::value,void>::type  
       real_factorial(unsigned int const n, T &xfact)
  { 
   T res;

   bool zdebug = false;

   std::ostringstream os;

   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "          >>>>> factorial(unsigned int)"
         << "\n\n"
         << "          n = " << n
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   //

   xfact = 0.0e+00;

   //

   if(0 == n)
     {
      res = 1.0e+00;
     }
   else if(n == 1)
     {
      res = 1.0e+00;
     }
   else if(n > 1)
     {
      int const n1 = n - 1;

      T ttemp;

      real_factorial(n1,ttemp);

      res = ( static_cast<T>(n) ) * ttemp;
     }
   else
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "          **** Error: real_factorial(unsigned int)"
         << "\n\n"
         << "          n is negative;  value = " << n
         << "\n\n";

      std::cout << os.str() << "\n";
     
      exit(0);
     }

   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "          Result: " << res 
         << "\n\n" 
         << "          <<<<< Completed: real_factorial(unsigned int)"
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   //

   xfact = res;

   return;  
  }
   // End function factorial()

/**
 *  Function: complex_unnormalized_assoc_regular_legendre()
 *
 *  This is the pictorial view of cpmvec[][]. It is a lower 
 *  half triangle of values. We build it using recursion 
 *  relationships along the diagonal firs and down each column
 *
 *           <---- m index ----->
 *               012345.....
 *             +
 *           0 |\
 *           1 | \
 *  l-index  2 |  \
 *           3 |   \
 *           4 |    \
 *           5 |     \
 *          ...|     |
 *           9 |     |
 *         
 *
 */
template <typename T> 
   typename std::enable_if<std::is_floating_point<T>::value,void>::type  
     complex_unnormalized_assoc_regular_legendre(
           unsigned int const mmax, 
           unsigned int const lmax,
           std::complex<T> z,
           std::vector<std::vector<std::complex<T> > > &cpmvec)
  {
   bool const zdebug = true;

   std::string const method_name_str = "complex_unnormalized_assoc_regular_legendre()";

   //
   //---- Debug banner header 
   //

   if(zdebug)
     {
      printf("\n\n");
      printf("      ====> %s <==== ", method_name_str.c_str());
      printf("\n\n");
      printf("      Input data: \n");
      printf("        Lmax (m) = %4u \n", lmax);
      printf("        Mmax (n) = %4u   ", mmax);
      printf("\n\n");

      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float)  ) )
        {
         cformat_str = "        z = (%10.6f,%10.6f) \n";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "        z = (%10.6Lf,%10.6Lf) \n";
        }

      printf(cformat_str.c_str(), z.real(), z.imag());

      printf("\n\n");
      printf("      **** End of input data");
      printf("\n\n");
     }

   //
   //---- Let's establish the constants 0.0, 1.0, 2.0 
   //     as complex variables
   //
   //     Following initialisation works over different versions
   //     of the compiler 
   //

   std::complex<T> ZERO; std::complex<T> ZONE; std::complex<T> ZTWO;

   ZERO.real(0.0e+00);   ZONE.real(1.0e+00);   ZTWO.real(2.0e+00);
   ZERO.imag(0.0e+00);   ZONE.imag(0.0e+00);   ZTWO.imag(0.0e+00);

   //
   //---- Initialise the output vectors 
   //

   for(int lval=0; lval<=lmax; ++lval)
      {
       for(int mval = 0; mval<=lval; ++mval)
          {
           cpmvec[lval][mval].real(0.0e+00);
           cpmvec[lval][mval].imag(0.0e+00);
          } 
      }

   //
   //---- Set the well known value
   //
   //         P_{00} (z) = 1.0
   //

   cpmvec[0][0] = ZONE;

   //
   //--- Are we at the point +1 or -1 on the real axis.
   //

   T const drealabs = z.real();
   T const dimagabs = z.imag();

   bool const z_is_real_and_abs_is_1 = (drealabs == 1.0e+00) && 
                                       (dimagabs == 0.0e+00);

   //

   T const zabs = std::abs(z);

   bool const z_inside_unit_circle  = zabs < 1.0e+00;

   bool const z_outside_unit_circle = zabs >= 1.0e+00;

   //

   if(zdebug)
     {
      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float)  ) )
        {
         cformat_str = "      Absolute value of z = %13.6f ";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "      Absolute value of z = %13.6Lf ";
        }

      printf(cformat_str.c_str(), zabs);
      printf("\n\n");

      if(z_is_real_and_abs_is_1)
        {
         printf("      Argument, z, is either +1 or -1");
         printf("\n\n");
        }

      if(z_outside_unit_circle)
        {
         printf("      Argument, z, lies on or outside the unit circle");
         printf("\n\n");
        }

      if(z_inside_unit_circle)
        {
         printf("      Argument, z, lies inside the unit circle");
         printf("\n\n");
        }
     }

   //
   //---- Cases:
   //
   //      1.   z = +1 or -1 therefore on real axis
   //
   //      2.   z is inside the unit circle
   //
   //      3.   z is on the unit circle or outside it
   //

   if(z_is_real_and_abs_is_1)
     {
      /**
      do i = 1, n
         cpm(0,i) = x ** i
         cpd(0,i) = 0.5D+00 * i * ( i + 1 ) * x ** ( i + 1 )
      end do

      do j = 1, n
         do i = 1, m
            if ( i == 1 ) then
               cpd(i,j) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
            else if ( i == 2 ) then
               cpd(i,j) = -0.25D+00 &
                 * ( j + 2 ) * ( j + 1 ) * j * ( j - 1 ) * x ** ( j + 1 )
             end if
         end do
      end do
      */

      //goto return_point; 
     }

   //
   //---- Inside or outside the unit circle, we use the
   //     same formula, but where a constant factor 
   //     is defined differently for these two cases.
   //
   //     Let's now set that factor.
   //

   std::complex<T> zls_factor; 

   if(z_inside_unit_circle)
     {
      zls_factor.real(1.0e+00);
      zls_factor.imag(0.0e+00);
     }
   else
     {
      zls_factor.real(-1.0e+00);
      zls_factor.imag(0.0e+00);
     }

   //
   //---- Next compute some auxiliary variables 
   //

   std::complex<T> const zsqd = z * z;

   std::complex<T> const zone_minus_zsqd = ZONE - zsqd;

   std::complex<T> const zs = zls_factor * zone_minus_zsqd;

   std::complex<T> const zq = std::sqrt(zs);  

   //
   //---- Apply the recursion formula along the diagnonal
   //
   //     We know P_{00} to start and we build P{11}, P{22}, ...
   //     in that order.
   //
   //     Note that we use "s" as a shorthand for zls_factor below
   //
   //       P_{ll}(z) = -s(2l - 1) sqrt( s (1 - z*z) ) P_{l-1,l-1}(z) 
   //
   //
   //    

   for(int ll = 1; ll<= mmax; ++ll)
      {
       T const dt = - ( 2.0e+00 * static_cast<double>(ll) - 1.0e+00 );
       
       std::complex<T> ztemp1;

       ztemp1.real(dt); ztemp1.imag(0.0e+00); 

       //

       std::complex<T> const ztemp2 = zls_factor * ztemp1; 

       std::complex<T> const ztemp3 = ztemp2 * zq;

       cpmvec[ll][ll] = ztemp3 * cpmvec[ll-1][ll-1];
      }
       // End of loop generating P_{ll} values. 

   //
   //---- With P{mm} known we can recurse up in l value 
   //
   //        P_{m+1,m}(x) = (2m + 1)xP_{m,m)(x)
   //

   for(int mval=0; mval<=mmax; ++mval)
      {
       if( mval < lmax )  // Dont exceed Lmax - no storage
         {
          std::complex<T> ztemp1 = z * cpmvec[mval][mval];

          T const dtemp2 = 2.0e+00 * static_cast<T>(mval) + 1.0e+00;

          std::complex<T> ztemp2;

          ztemp2.real(dtemp2); ztemp2.imag(0.0e+00); 

          cpmvec[mval+1][mval] = ztemp2 * ztemp1;
         }
      }

   //
   //---- Now fill in all other P_{l,m} values 
   //
   //     This uses equation (4.4.11) in the book 
   //
   //     We deal with the terms inside the brackets []
   //     and apply the factor at the end.
   //
 
   for(int mval=0; mval<mmax; ++mval)
      {
       for(int lval=mval+2; lval<=lmax; ++lval)
          {
           //
           //.... First part of recursion 
           //
           //       (2*l - 1) z P_{l-1,m} (z) 
           //

           std::complex<T> zpart1;

           {
            T const dtemp11 = 2.0e+00 * static_cast<double>(lval) - 1.0e+00;

            std::complex<T> ztemp11;

            ztemp11.real(dtemp11); ztemp11.imag(0.0e+00);

            //

            std::complex<T> ztemp12 = z * cpmvec[lval-1][mval];

            zpart1 = ztemp11 * ztemp12;
           }
            // End scope 1

           //
           //.... Second part of recursion
           //
           //       (l + m - 1) P_{l-2,m} 
           //

           std::complex<T> zpart2;           

           {
            T const dtemp21 = static_cast<T>(lval) + 
                              static_cast<T>(mval) - 1.0e+00;

            std::complex<T> ztemp21;

            ztemp21.real(dtemp21); ztemp21.imag(0.0e+00);

            zpart2 = ztemp21 * cpmvec[lval-2][mval];
           }
            // End scope 2

           //
           //.... Prepare the term in brackets in equation (4.4.11)
           //

           std::complex<T> zbrackets = zpart1 - zpart2;

           //

           std::complex<T> zfactor;

           {
            int const itemp = lval - mval;

            T const dtemp3 = static_cast<T>(itemp);  

            T const dtemp4 = 1.0e+00/dtemp3;

            zfactor.real(dtemp4); zfactor.imag(0.0e+00);
           }
            // End scope for building zfactor

           //
           //.... Compete the computation for equation (4.4.11)
           //

           cpmvec[lval][mval] = zfactor * zbrackets; 
          }
           // End loop over lval
      }
       // End loop over mval

   //
   //---- Return point 
   //

   if(zdebug)
     {
      printf("\n\n");
      printf("      **** Completed - %s ", method_name_str.c_str());
      printf("\n\n");
     }

   return;
  }
   // End of complex_unnormalized_assoc_regular_legendre()

/**
 *  Function: complex_unnormalized_assoc_irregular_legendre()
 *
 */
template <typename T> 
   typename std::enable_if<std::is_floating_point<T>::value,void>::type  
     complex_unnormalized_assoc_irregular_legendre(
           unsigned int const mmax, 
           unsigned int const lmax,
           std::complex<T> z,
           std::vector<std::vector<std::complex<T> > > &zcqmvec)
  {
   bool const zdebug = true;

   std::string const method_name_str = "complex_unnormalized_assoc_irregular_legendre()";

   //
   //---- Debug banner header 
   //

   if(zdebug)
     {
      printf("\n\n");
      printf("      ====> %s <==== ", method_name_str.c_str());
      printf("\n\n");
      printf("      Input data: \n");
      printf("        Lmax (n) = %4u \n", lmax);
      printf("        Mmax (m) = %4u   ", mmax);
      printf("\n\n");

      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float)  ) )
        {
         cformat_str = "        z = (%10.6f,%10.6f) \n";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "        z = (%10.6Lf,%10.6Lf) \n";
        }

      printf(cformat_str.c_str(), z.real(), z.imag());

      printf("\n\n");
      printf("      **** End of input data");
      printf("\n\n");
     }

   //
   //---- Let's establish the constants 0.0, 1.0, 2.0 
   //     as complex variables
   //
   //     Following initialisation works over different versions
   //     of the compiler 
   //
  
   std::complex<T> ZMINUS1;

   ZMINUS1.real(-1.0e+00);
   ZMINUS1.imag(0.0e+00);

   std::complex<T> ZERO;  std::complex<T> ZONE;  std::complex<T> ZTWO;

   ZERO.real(0.0e+00);         ZONE.real(1.0e+00);         ZTWO.real(2.0e+00);
   ZERO.imag(0.0e+00);         ZONE.imag(0.0e+00);         ZTWO.imag(0.0e+00);

   //

   std::complex<T> ZHALF;

   ZHALF.real(0.5e+00); 
   ZHALF.imag(0.0e+00);

   //
   //---- Initialise the output vectors 
   //

   for(int lval=0; lval<=lmax; ++lval)
      {
       for(int mval = 0; mval<=mmax; ++mval)
          {
           zcqmvec[lval][mval].real(0.0e+00);
           zcqmvec[lval][mval].imag(0.0e+00);
          } 
      }

   //
   //--- Are we at the point +1 or -1 on the real axis.
   //

   T const drealabs = z.real();
   T const dimagabs = z.imag();

   bool const z_is_real_and_abs_is_1 = (drealabs == 1.0D+00) && 
                                       (dimagabs == 0.0D+00);

   //

   T const zabs = std::abs(z);

   bool const z_inside_unit_circle  = zabs < 1.0e+00;

   bool const z_outside_unit_circle = zabs >= 1.0e+00;

   //

   if(zdebug)
     {
      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float) ) )
        {
         cformat_str = "      Absolute value of z = %13.6f ";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "      Absolute value of z = %13.6Lf ";
        }

      printf(cformat_str.c_str(), zabs);
      printf("\n\n");

      //

      if(z_is_real_and_abs_is_1)
        {
         printf("      Argument, z, is either +1 or -1");
         printf("\n\n");

         return;
        }

      if(z_outside_unit_circle)
        {
         printf("      Argument, z, lies on or outside the unit circle");
         printf("\n\n");
        }

      if(z_inside_unit_circle)
        {
         printf("      Argument, z, lies inside the unit circle");
         printf("\n\n");
        }
     }

   //
   //---- Inside or outside the unit circle, we use the
   //     same formula, but where a constant factor 
   //     is defined differently for these two cases.
   //
   //     Let's now set that factor.
   //

   std::complex<T> zls_factor; 

   if(z_inside_unit_circle)
     {
      zls_factor.real(1.0e+00);
      zls_factor.imag(0.0e+00);
     }
   else
     {
      zls_factor.real(-1.0e+00);
      zls_factor.imag(0.0e+00);
     }

   //
   //---- Next compute some auxiliary variables 
   //

   std::complex<T> const zsqd = z * z;

   std::complex<T> const zone_minus_zsqd = ZONE - zsqd;

   std::complex<T> const zs = zls_factor * zone_minus_zsqd;

   std::complex<T> const zq = std::sqrt(zs);  

   //
   //---- Computation of zcq0 
   //
   //      (1/2) * log (s * (z + 1) / (z - 1))
   //
   //     The factor s represent zls_factor above.
   //

   std::complex<T> zcq0; 

   {
    std::complex<T> one_plus_z  = ZONE + z; 

    std::complex<T> one_minus_z = ZONE - z;

    std::complex<T> ztemp1 = one_plus_z / one_minus_z;

    std::complex<T> ztemp2 = zls_factor * ztemp1;

    std::complex<T> ztemp3 = std::log(ztemp2); 

    zcq0 = ZHALF * ztemp3;
   }
    // End of scope for computation of zcq0

   if(zdebug)
     {
      T testprint = 0.0e+00;

      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float) ) )
        {
         cformat_str = "      zcq0 = (%13.6e,%13.6e) "; 
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "      zcq0 = (%13.6Le,%13.6Le) "; 
        }

      printf(cformat_str.c_str(), zcq0.real(), zcq0.imag()); 
      printf("\n\n");
     }

/**
  if ( abs ( x ) == 1.0D+00 .and. y == 0.0D+00 ) then
    do i = 0, m
      do j = 0, n
        cqm(i,j) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
        cqd(i,j) = cmplx ( 1.0D+30, 0.0D+00, kind = 8 )
      end do
    end do
    return
  end if

  cq0 = 0.5D+00 * log ( ls * ( 1.0D+00 + z ) / ( 1.0D+00 - z ) )

  if ( xc < 1.0001D+00 ) then

    cqm(0,0) = cq0
    cqm(0,1) = z * cq0 - 1.0D+00
    cqm(1,0) = -1.0D+00 / zq
    cqm(1,1) = - zq * ( cq0 + z / ( 1.0D+00 - z * z ) )
    do i = 0, 1
      do j = 2, n
        cqm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * z * cqm(i,j-1) &
          - ( j + i - 1.0D+00 ) * cqm(i,j-2) ) / ( j - i )
      end do
    end do

    do j = 0, n
      do i = 2, m
        cqm(i,j) = -2.0D+00 * ( i - 1.0D+00 ) * z / zq * cqm(i-1,j) &
          - ls * ( j + i - 1.0D+00 ) * ( j - i + 2.0D+00 ) * cqm(i-2,j)
      end do
    end do

  else
*/
   
   //
   //---- 
   //

   if(z_outside_unit_circle)
     {
      unsigned int km;

      if(1.1D+00 < zabs)
        {
         km = 40 + lmax + mmax;
        }
      else
        {
         double const dtemp1 = std::log(zabs - 1.0e+00);

         double const dtemp2 = -1.0e+00 - 1.8e+00*dtemp1; 

         km = 40 + lmax + mmax;

         km *=  static_cast<double>(floor(dtemp2)); // NB: floor type double arg
        }

      if(zdebug)
        {
         printf("      For the backwards recursion, km = %d ", km);
         printf("\n\n");
        }

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //
      //.... Iterate backwards to find F_{k}^{0}
      //
      //     This means that we apply equation (4.5.10) 
      //     in the book:
      //
      //                   (2k+3)zF_{k+1}^{0} - (k+2)F_{k+2}^{0}
      //       F_{k}^{0} = ------------------------------------- 
      //                               (k+1)
      //
      //     At the end of the next for loop, we have set
      //
      //       Q_{k}^{0} (z) = F_{k}^{0} (z) 
      //
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      std::complex<T> zcqf2 = ZERO;
      std::complex<T> zcqf1 = ZONE;
      std::complex<T> zcqf0 = ZERO;

      for(int k = km; k>= 0; k--)
         {
          //
          //.... Ok, build the left hand side term within the bracket
          //

          int const itwok = 2 * k;

          int const itwok_plus_three = itwok + 3;

          T const dtwok_plus_three = static_cast<T>(itwok_plus_three);

          std::complex<T> ztwok_plus_three;

          ztwok_plus_three.real(dtwok_plus_three); 
          ztwok_plus_three.imag(0.0e+00);

          std::complex<T> ztemp1 = ztwok_plus_three * z;

          std::complex<T> znumerator1 = ztemp1 * zcqf1; 

          //
          //.... Ok, build the right hand side term within the bracket
          //

          int const k_plus_two     = k + 2;

          T const dk_plus_two = static_cast<T>(k_plus_two);

          std::complex<T> zkplus_two;

          zkplus_two.real(dk_plus_two);
          zkplus_two.imag(0.0e+00);

          std::complex<T> znumerator2 = zkplus_two * zcqf2; 

          //
          //.... Compute the bracket in equation (4.5.10)
          //
          //     Its the numerator of the fraction
          //

          std::complex<T> znumerator = znumerator1 - znumerator2;

          //
          //.... Compute the denominator of equation (4.5.10)
          //

          int const k_plus_one = k + 1;

          T const dk_plus_one = static_cast<T>(k_plus_one);

          std::complex<T> zk_plus_one;

          zk_plus_one.real(dk_plus_one);
          zk_plus_one.imag(0.0e+00);

          std::complex<T> zdenominator = zk_plus_one;

          //
          //.... All the parts are in place
          //
          //     Compute F_{k}^{0} in equation (4.5.10)
          //
          //     We only store the F() values whenever
          //     k <= lmax 
          //

          zcqf0 = znumerator / zdenominator;

          if( k <= lmax ) 
            {
             zcqmvec[k][0] = zcqf0;
            }

          //
          //.... Get ready for the next iteration by "pushing"
          //     the values for the "next" k.
          //

          zcqf2 = zcqf1;
          zcqf1 = zcqf0;
         }
          // End of loop over k - equation (4.5.10)

      //
      //---- For all "l" values we multiply by a factor, so that
      //     at the end of this next loop we have 
      //
      //                                        Q_{0}^{0} (z)
      //       Q_{k}^{0} (z) = F_{k}^{0} (z) * ---------------
      //                                        F_{0}^{0} (z)  
      //

      std::complex<T> zfactor = zcq0 / zcqf0; 

      for(int k=0; k<=lmax; ++k)
         {
          std::complex<T> ztemp = zcqmvec[k][0] * zfactor; 

          zcqmvec[k][0] = ztemp; 
         }

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //
      //.... Iterate backwards to find F_{k}^{1}
      //
      //     This means that we apply equation (4.5.10) 
      //     in the book, as above but this time to find 
      //
      //                       (2k+3)zF_{k+1}^{1} - (k+1)F_{k+2}^{1}
      //       F_{k}^{1} (z) = ------------------------------------- 
      //                                    (k+2)
      //
      //     At the end of this next for loop we have set:
      //     
      //       Q_{k}^{1} (z) = F_{k}^{1} (z) 
      //
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //
      //---- Reset the push down stack variables 
      //

      zcqf2 = ZERO;
      zcqf1 = ZONE;
      zcqf0 = ZERO;

      for(int k = km; k>= 0; k--)
         {
          //
          //.... Ok, build the lefthand side term within the bracket
          //

          int const itwok = 2 * k;

          int const itwok_plus_three = itwok + 3;

          T const dtwok_plus_three = static_cast<T>(itwok_plus_three);

          std::complex<T> ztwok_plus_three;

          ztwok_plus_three.real(dtwok_plus_three); 
          ztwok_plus_three.imag(0.0e+00);

          std::complex<T> ztemp1 = ztwok_plus_three * z;

          std::complex<T> znumerator1 = ztemp1 * zcqf1; 

          //
          //.... Ok, build the righthand side term within the bracket
          //


          int const k_plus_one     = k + 1;

          T const dk_plus_one = static_cast<T>(k_plus_one);

          std::complex<T> zkplus_one;

          zkplus_one.real(dk_plus_one);
          zkplus_one.imag(0.0e+00);

          std::complex<T> znumerator2 = zkplus_one * zcqf2; 

          //
          //.... Compute the bracket in equation (4.5.10)
          //
          //     Its the numerator of the fraction
          //

          std::complex<T> znumerator = znumerator1 - znumerator2;

          //
          //.... Compute the denominator of equation (4.5.10)
          //

          int const k_plus_two = k + 2;

          T const dk_plus_two = static_cast<T>(k_plus_two);

          std::complex<T> zk_plus_two;

          zk_plus_two.real(dk_plus_two);
          zk_plus_two.imag(0.0e+00);

          std::complex<T> zdenominator = zk_plus_two;

          //
          //.... All the parts are in place
          //
          //     Compute F_{k}^{0} in equation (4.5.10)
          //
          //     We only store the F() values whenever
          //     k <= lmax 
          //

          zcqf0 = znumerator / zdenominator;

          if( k <= lmax ) 
            {
             zcqmvec[k][1] = zcqf0;
            }

          //
          //.... Get ready for the next iteration by "pushing"
          //     the values for the "next" k.
          //

          zcqf2 = zcqf1;
          zcqf1 = zcqf0;
         }
          // End loop over k using equation (4.5.10)  

      //
      //---- For all "l" values we multiply by a factor, so that
      //     at the end of this next loop we have 
      //
      //                                                -1
      //     Q_{k}^{1} (z) = F_{k}^{1} (z) * --------------------------
      //                                     sqrt(s(1-z*z))*F_{0}^{1}(z)  
      //

      std::complex<T> zfactor1_temp = ZMINUS1 / zq;
 
      std::complex<T> zfactor1 = zfactor1_temp / zcqf0; 

      for(int k=0; k<=lmax; ++k)
         {
          std::complex<T> ztemp = zcqmvec[k][1] * zfactor1; 

          zcqmvec[k][1] = ztemp; 
         }
      
      //
      //.... Given that we have computed Q_{l}_{0} and Q_{l}^{m}
      //     for all l values up to lmax, we now apply the
      //     equation (4.4.10) to generate all "m" values up to
      //     mmax. The recursion is show below:
      //
      //                  -2 (m-1)
      //     Q_{l}^{m} = ---------- * z * Q_{l}^{m-1} -  
      //                 sqrt(1-z*z)
      //
      //
      //                   (l + m - 1)(l - m + 2) * Q_{l}^{m-2}
      //

      for(int lval=0; lval<=lmax; ++lval)
         {
          std::complex<T> zcq0 = zcqmvec[lval][0];
          std::complex<T> zcq1 = zcqmvec[lval][1];

          for(int mval=2; mval<=mmax; ++mval)
             {
              //
              //.... First term in the recursion of (4.4.10)
              //

              std::complex<T> zterm_first = ZERO;

              {
               T const dmval_plus_1 = static_cast<T>(mval - 1);

               T const dtemp11      = -2.0e+00 * dmval_plus_1;  

               std::complex<T> ztemp11;

               ztemp11.real(dtemp11); ztemp11.imag(0.0e+00);

               std::complex<T> ztemp12 = ztemp11 * z;

               std::complex<T> ztemp13 = ztemp12 / zq;

               zterm_first  = ztemp13 * zcq1;
              }
               // End of scoping for computation of first term in (4.4.10)

              //
              //.... Second term in the recursion of (4.4.10)
              //

              std::complex<T> zterm_second = ZERO;

              {
               T const dlval_plus_mval_minus_1 
                              = static_cast<T>(lval + mval + 1);

               T const dlval_minus_mval_plus_2
                              = static_cast<T>(lval - mval + 2);

               T const dtemp21 = dlval_plus_mval_minus_1 * 
                                        dlval_minus_mval_plus_2;

               std::complex<T> ztemp21;

               ztemp21.real(dtemp21); ztemp21.imag(0.0e+00);

               zterm_second  = ztemp21 * zcq0;
              }
               // End of scoping for computation of second term in (4.4.10)

              //
              //.... Compute the recursion and store 
              //
              //     Note that its a subtraction
              //

              std::complex<T> zcqf = zterm_first - zterm_second;

              zcqmvec[lval][mval] = zcqf;

              //
              //.... Push down the values for the next term
              //

              zcq0 = zcq1;
              zcq1 = zcqf;
             }
              // End loop over m values 
         }
          // End loop over l values
     }
      // End of computation for z lying outside the unit circle

   //
   //---- Return point 
   //

   if(zdebug)
     {
      printf("\n\n");
      printf("      **** Completed - %s ", method_name_str.c_str());
      printf("\n\n");
     }

   return;
  }
   // End of ...()

/**
 *  Method: power_complx_to_unsigned_int()
 *
 *  Given a complex number, raise it to a positive integer power. 
 *
 *  Restrict instantiations to std::complex<float>, 
 *                             std::complex<double>      and
 *                             std::complex<long double>
 *
 */
template <typename T, bool ZDEBUG> 
   std::complex< typename std::enable_if<std::is_floating_point<T>::value,T>::type > 
           power_complx_to_unsigned_int(std::complex<T> const &z, 
                                        unsigned int    const ipow)
  {
   if(ZDEBUG)
     {
      std::ostringstream os;

      os.str(""); os.clear();

      os << "\n\n"
         << "          >>>>> Template instantiation: power_complx_to_unsigned_int()"
         << "\n\n"
         << "          z = (" << z.real() << ", " << z.imag() << ") " << " - power (ipow) = " << ipow
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   //

   std::complex<T> const zsqd = z*z;

   std::complex<T> zpow(0.0e+00,0.0e+00);
   
   //

   switch(ipow)
     {
      case  0: 
              zpow.real(1.0e+00);
              zpow.imag(0.0e+00);
              break;

      case  1:
              zpow = z;
              break;

      case 2: 
              zpow = z*z;
              break;

      case 3: 
              zpow = z*z*z;
              break;

      case 4: 
              zpow = zsqd * zsqd;
              break;

      default:
              zpow = z;

              for(int i=2; i<=ipow; i++)
                 {
                  zpow = zpow * z;
                 }
              break;
     }

   //
   //---- Return point
   //

   if(ZDEBUG)
     {
      std::ostringstream os;

      os.str(""); os.clear();

      os << "          Result: zpow = (" << zpow.real() << ", " << zpow.imag() << ") " 
         << "\n\n" 
         << "          <<<<< Completed: power_complx_to_unsigned_int()"
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   return zpow;
  }
   // End of template power_complx_to_unsigned_int()

/**
 *    cplm()
 *                                                                       
 *    Calculates the regular associated legendre function  P  ( z )    
 *                                                          l,m       
 *    l,m assumed to be integer                            
 *    z   complex and |z| > 1.0 
 *                    
 *    Developed from an initial version by Bell and McLaughlin (1983)
 *    Never published.          
 *                                      
 */
template <typename T, bool ZDEBUG> 
   typename std::enable_if<std::is_floating_point<T>::value,void>::type  
     cplm(unsigned int const l, 
          unsigned int const m,
          std::complex<T> const z,
          std::complex<T> zretVal)                                  
  {
   //
   //.... Local complex constants  
   //
   
   std::complex<double> const ZERO(0.0e+00, 0.0e+00); 
   std::complex<double> const ZONE(1.0e+00, 0.0e+00); 
   std::complex<double> const ZTWO(2.0e+00, 0.0e+00);

   //
   //.... Local double precision constants 
   //

   double const  DTHRESHOLD_POLE = 1.0e-09; 

   double const  DHALF = 0.5d+00;
   double const  DONE  = 1.0d+00;
   double const  DTWO  = 2.0d+00;

   //
   //.... Local boolean variables  
   //

   bool const zdebug = true;

   //
   //....  Output string stream 
   //

   std::ostringstream os;

   //
   //---- Debug banner header 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "     >>>> cplm - compute associated Legendre functions "
         << "\n\n"
         << "     Input data: "  
         << " l = " << l 
         << " m = " << m 
         << " z = " << z 
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   //
   //---- Default the return value
   //

   zretVal = ZERO;

   //
   //---- Check arguments
   //
   //     l should be >= 0  
   //     m should be [-l,+l]
   //
   //     z should not be (1,0) - need to pick a delta neighbouroud.
   //
   //   if(l .lt. 0)then
   //     write(iwrite,9000)
   //     stop 999
   //   endif
   //
   //   if((m .lt. -l) .or. (m > l))then
   //     write(iwrite,9000)
   //     stop 999
   //   endif 
   //
   //.... Distance from (1,0) - the pole 
   //

   std::complex<double> zdiff = z - ZONE;

   double dt = std::abs(z);

   if(dt < DTHRESHOLD_POLE)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     ***** Error in cplm - compute associated Legendre functions"
         << "\n\n"
         << "     Z is too close to the singularity at +1 on real axis " 
         << "\n\n";

      std::cout << os.str() << "\n";

      exit(-1);
     } 

   //
   //---- Prepare values for use in computation
   //

   int    lm  = l - m;                       

   int    lm2 = lm/2;                                                    

   //
   //---- First deal with term (ir = 0)
   //
   //               (2 * l)!
   //        ---------------------------   *  z^(l - m)
   //          (2^l) * (l!) * (l - m)! 
   //
   //     The pow() function gets bad press seeingly for int
   //     arguments, since it works with real numbers. So, let's
   //     roll our own here.
   //

   int const  itwo_l = 2*l;

   //

   if(l < 0)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     ***** Error in cplm " 
         << "\n\n"
         << "     \"l\" is negative. "  
         << "\n\n";

      std::cout << os.str() << "\n";

      exit(0);
     }

   //

   int itwo_power_l = 0;

   if(0 == l)
     {
      itwo_power_l = 1;
     }
   else if(1 == l)
     {
      itwo_power_l = 2;
     }
   else if( (l>=2) && (l<=31) )
     {
      itwo_power_l = 2 << (l - 1);
     }
   else
     {
      exit(1);
     }

   double const xtwo_power_l = static_cast<double>(itwo_power_l);

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     l = " << l << ",  2^l = " << itwo_power_l << " and as float = " << xtwo_power_l
         << "\n\n";

      std::cout << os.str() << "\n";

     } 

   //

   double dfactorial_l, dfactorial_lm, dfactorial_2l;

   real_factorial(l,dfactorial_l);

   real_factorial(lm,dfactorial_lm);

   real_factorial(itwo_l,dfactorial_2l);

   //

   double const denominator = xtwo_power_l * dfactorial_l * dfactorial_lm;                  

   double       wlmro       = dfactorial_2l / denominator;                  

   std::complex<double> zwlmro;

   zwlmro.real(wlmro);  zwlmro.imag(0.0e+00);

   //

   std::complex<double> zpow = power_complx_to_unsigned_int<double,true>(z,lm);

   std::complex<double> zsum = zwlmro * zpow;

   //
   //---- Now deal with terms ir = 1, 2, ... lm2 - indeed if any 
   //                                  
   //     We compute the following term, "wlmr",  at each value "ir" 
   //
   //            -1 * (l - m - 2*ir + 2) * (l - m + 1)
   //           --------------------------------------- * z^(l - m - 2*ir)
   //                 (2*ir) * 2 * (l - ir + 1)
   //
   //     Note the product relationship as we save this to "wlmro" and reuse
   //     at next iteration.
   //
 
   for(int ir=1; ir<=lm2; ++ir)
      {
       int const lmir = lm - 2*ir;                                                 
       int const jlir = l  - ir;                                                  

       double wlmr = -wlmro*(static_cast<double>(lmir) + DTWO)*(static_cast<double>(lmir) + DONE);                               
       
       wlmr = wlmr/(DTWO*static_cast<double>(ir)*(DTWO*static_cast<double>(jlir)+ DONE));                

       //

       std::complex<double> zwlmr;

       zwlmr.real(wlmr); zwlmr.imag(0.0e+00);

       //

       std::complex<double> const zpowterm = power_complx_to_unsigned_int<double,true>(z,lmir);

       std::complex<double> const t = zwlmr * zpowterm;                                             

       zsum  = zsum + t;                                                    
         
       wlmro = wlmr;                                                  
      }                                                         

   //
   //---- Now set the final value into the variable "p"
   //

   if(m == 0)
     {                                             
      zretVal = zsum;
     }                               
   else
     {                                             
      double const xn  = static_cast<double>(m) * DHALF;                                          

      std::complex<double> const zsqd   = z*z;

      std::complex<double> const ztemp  = zsqd - ZONE;

      std::complex<double> const ztemp2 = std::pow(ztemp, xn);

      zretVal = ztemp2 * zsum;                             
     }  

   //
   //---- Return point 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "     Computed value of associated Legendre polynomial = " 
         << "(" << zretVal.real() << ", " << zretVal.imag() << ") " 
         << "\n\n"
         << "     ***** cplm - completed \n";

      std::cout << os.str() << "\n";
     }

   return;
  }
   // End function cplm() 

#endif // For _ASSOC_LEGENDRE_TEMPLATED_COMPLEX_INCLUDE_176496_H_

//********************************************************************** 
//********************************************************************** 
//
//   End of file 
//
//********************************************************************** 
//********************************************************************** 
