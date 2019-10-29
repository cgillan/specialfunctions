//************************************************************************
//************************************************************************
//
//   File: zhang_jin.cxx  
//
//   This file implements the equations for the computation of
//   associated Legendre polynomials of the first and second 
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

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <cmath>
#include <cfenv>
#include <climits> 
#include <cfloat> 
#include <vector>
#include <algorithm>
#include <complex>
#include <cstdint>

//#include "associated_legendre_functions.hxx"
  
void complex_unnormalized_assoc_regular_legendre(
           unsigned int const m, 
           unsigned int const n,
           std::complex<double> z,
           std::vector<std::vector<std::complex<double> > > &cpmvec);

void complex_unnormalized_assoc_irregular_legendre(
           unsigned int const mmax, 
           unsigned int const lmax,
           std::complex<double> z,
           std::vector<std::vector<std::complex<double> > > &zcqmvec);

#pragma STDC FENV_ACCESS ON

/**
 *   Main program - test harness
 *
 */

int main(int argc, char **argv)
  {
   int const Lmax = 5;

   int const Mmax = 5;

   std::complex<double> z;

   //===================================================================
   //
   //   T E S T   R E G U L A R   A S S O C I A T E D   L E G E N D R E 
   //
   //===================================================================

   std::cout << "      Numeric limits   \n";
   std::cout << "      -------------- \n\n";
   std::cout << "      int:             \n"; 
   std::cout << "        Minimum value: " << std::numeric_limits<int>::min()        << '\n';
   std::cout << "        Maximum value: " << std::numeric_limits<int>::max()        << '\n';
   std::cout << "        int is signed: " << std::numeric_limits<int>::is_signed    << '\n';
   std::cout << "        Non-sign bits: " << std::numeric_limits<int>::digits       << '\n';
   std::cout << "        has infinity:  " << std::numeric_limits<int>::has_infinity << '\n';
   std::cout << "      float:           \n";
   std::cout << "        minimum value : " << FLT_MIN        << "\n";
   std::cout << "        maximum value : " << FLT_MAX        << "\n";
   std::cout << "        epsilon       : " << FLT_EPSILON    << "\n";
   std::cout << "        exponent radix: " << FLT_RADIX      << "\n";
   std::cout << "      double:           \n";
   std::cout << "        minimum value : " << DBL_MIN        << "\n";
   std::cout << "        maximum value : " << DBL_MAX        << "\n";
   std::cout << "        epsilon       : " << DBL_EPSILON    << "\n";
   std::cout << "      long double:      \n";
   std::cout << "        minimum value : " << LDBL_MIN       << "\n";
   std::cout << "        maximum value : " << LDBL_MAX       << "\n";
   std::cout << "        epsilon       : " << LDBL_EPSILON   << "\n";

   //===================================================================
   //
   //   T E S T   R E G U L A R   A S S O C I A T E D   L E G E N D R E 
   //
   //===================================================================

   std::vector<std::vector<std::complex<double> > > cpmvec;

   cpmvec.resize(Lmax+1);

   for(int lval=0; lval<=Lmax; ++lval)
      {
       cpmvec[lval].resize(lval+1);

       for(int mval=0; mval<=lval; ++mval)
          {
           cpmvec[lval][mval].real(0.0e+00);
           cpmvec[lval][mval].imag(0.0e+00);
          }
      }

   //

   z.real(2.5e+00); z.imag(0.0e+00);

   //
   //---- Clear all floating point exceptions
   //     before calling the compute method and 
   //     then test after the call for such 
   //     exceptions.
   //
 
   std::feclearexcept(FE_ALL_EXCEPT);

   complex_unnormalized_assoc_regular_legendre(Lmax, Lmax, z, cpmvec);

   if(std::fetestexcept(FE_ALL_EXCEPT))
     {
      if(std::fetestexcept(FE_DIVBYZERO)) 
        {
         printf("      **** Error: division by zero reported during computation    \n");
         printf("                  of the irregular associated Legendre functions.   ");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_INEXACT)) 
        {
         printf("      **** Error: inexact computation reported during computation \n");
         printf("                  of the irregular associated Legendre functions.   ");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_INVALID)) 
        {
         printf("      **** Error: invalid computation reported during computation \n");
         printf("                  of the irregular associated Legendre functions.   ");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_OVERFLOW)) 
        {
         printf("      **** Error: overflow reported during computation of the \n");
         printf("                  irregular associated Legendre functions.      ");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_UNDERFLOW)) 
        {
         printf("      **** Error: underflow reported during computation of the \n");
         printf("                  irregular associated Legendre functions.       ");
         printf("\n\n");
        }
     }
      // End of test on floating point exceptions being raised

   //

   printf("\n\n");
   printf("     Complex argument (z) = (%15.6f,%15.6f) ", z.real(), z.imag());
   printf("\n\n");
   printf("      Computed associated Legendre functions of the first kind (polynomials)");
   printf("\n\n");
   printf("      Note that these values are complex numbers in general");
   printf("\n\n");

   for(int m=0; m<=Lmax; ++m)
      {
       for(int l=m; l<=Lmax; ++l)
          {
           printf("      l=%2d m=%2d   (%13.6e, %13.6e) \n", l, m, cpmvec[l][m].real(), cpmvec[l][m].imag());
          }
      }

   //=====================================================================
   //
   //   T E S T  I R R E G U L A R   A S S O C I A T E D   L E G E N D R E 
   //
   //=====================================================================

   std::vector<std::vector<std::complex<double> > > cqmvec;

   cqmvec.resize(Lmax+1);

   for(int lval=0; lval<=Lmax; ++lval)
      {
       cqmvec[lval].resize(Mmax+1);

       for(int mval=0; mval<=Mmax; ++mval)
          {
           cqmvec[lval][mval].real(0.0e+00);
           cqmvec[lval][mval].imag(0.0e+00);
          }
      }

   //

   z.real(2.5e+00); z.imag(0.0e+00);

   //

   complex_unnormalized_assoc_irregular_legendre(Mmax, Lmax, z, cqmvec);

   //

   printf("\n\n ");

   for(int icol=2; icol <=72; ++icol) printf("-");

   printf("\n\n");
   printf("     Complex argument (z) = (%15.6f,%15.6f) ", z.real(), z.imag());
   printf("\n\n");
   printf("     Computed associated Legendre functions of the second kind (irregular)");
   printf("\n\n");
   printf("     Note that these values are complex numbers in general");
   printf("\n\n");

   printf("     index  l    m               z                        |z|        Associated Legendre function \n");
   printf("     ----- ---  ---  ------------------------------  -------------  ----------------------------- \n");

   for(int m=0; m<=Mmax; ++m)
      {
       for(int l=0; l<=Lmax; ++l)
          {
           printf("     %5d %3d  %3d  (%13.6e, %13.6e)  %13e  (%13e,%13e) \n",  
             0, l, m, 
             z.real(), z.imag(),  
             std::abs(z),
              cqmvec[l][m].real(), cqmvec[l][m].imag());
          }
      }

   printf("\n\n");
   printf("     **** Completed computations ");
   printf("\n\n");
  }
   // End of main

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
void complex_unnormalized_assoc_regular_legendre(
           unsigned int const m, 
           unsigned int const n,
           std::complex<double> z,
           std::vector<std::vector<std::complex<double> > > &cpmvec)
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
      printf("        Lmax (m) = %4u \n", m);
      printf("        Mmax (n) = %4u   ", n);
      printf("\n\n");
      printf("        z = (%10.6f,%10.6f) \n", z.real(), z.imag());
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

   std::complex<double> ZERO;  std::complex<double> ZONE;  std::complex<double> ZTWO;

   ZERO.real(0.0e+00);         ZONE.real(1.0e+00);         ZTWO.real(2.0e+00);
   ZERO.imag(0.0e+00);         ZONE.imag(0.0e+00);         ZTWO.imag(0.0e+00);

   //
   //---- Initialise the output vectors 
   //

   for(int lval=0; lval<= n; ++lval)
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

   double const drealabs = z.real();
   double const dimagabs = z.imag();

   bool const z_is_real_and_abs_is_1 = (drealabs == 1.0D+00) && 
                                       (dimagabs == 0.0D+00);

   //

   double const zabs = std::abs(z);

   bool const z_inside_unit_circle  = zabs < 1.0e+00;

   bool const z_outside_unit_circle = zabs >= 1.0e+00;

   //

   if(zdebug)
     {
      printf("      Absolute value of z = %13.6f ", zabs);
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

   std::complex<double> zls_factor; 

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

   std::complex<double> const zsqd = z * z;

   std::complex<double> const zone_minus_zsqd = ZONE - zsqd;

   std::complex<double> const zs = zls_factor * zone_minus_zsqd;

   std::complex<double> const zq = std::sqrt(zs);  

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

   for(int ll = 1; ll<= m; ++ll)
      {
       double const         dt = - ( 2.0e+00 * static_cast<double>(ll) - 1.0e+00 );
       
       std::complex<double> ztemp1;

       ztemp1.real(dt); ztemp1.imag(0.0e+00); 

       //

       std::complex<double> const ztemp2 = zls_factor * ztemp1; 

       std::complex<double> const ztemp3 = ztemp2 * zq;

       cpmvec[ll][ll] = ztemp3 * cpmvec[ll-1][ll-1];
      }
       // End of loop generating P_{ll} values. 

   //
   //---- With P{mm} known we can recurse up in l value 
   //
   //        P_{m+1,m}(x) = (2m + 1)xP_{m,m)(x)
   //

   for(int mval=0; mval<=m; ++mval)
      {
       if( mval < n )  // Dont exceed Lmax - no storage
         {
          std::complex<double> ztemp1 = z * cpmvec[mval][mval];

          double const dtemp2 = 2.0e+00 * static_cast<double>(mval) + 1.0e+00;

          std::complex<double> ztemp2;

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
 
   for(int mval=0; mval<m; ++mval)
      {
       for(int lval=mval+2; lval<=n; ++lval)
          {
           //
           //.... First part of recursion 
           //
           //       (2*l - 1) z P_{l-1,m} (z) 
           //

           std::complex<double> zpart1;

           {
            double const dtemp11 = 2.0e+00 * static_cast<double>(lval) - 1.0e+00;

            std::complex<double> ztemp11;

            ztemp11.real(dtemp11); ztemp11.imag(0.0e+00);

            //

            std::complex<double> ztemp12 = z * cpmvec[lval-1][mval];

            zpart1 = ztemp11 * ztemp12;
           }
            // End scope 1

           //
           //.... Second part of recursion
           //
           //       (l + m - 1) P_{l-2,m} 
           //

           std::complex<double> zpart2;           

           {
            double const dtemp21 = static_cast<double>(lval) + 
                                   static_cast<double>(mval) - 1.0e+00;

            std::complex<double> ztemp21;

            ztemp21.real(dtemp21); ztemp21.imag(0.0e+00);

            zpart2 = ztemp21 * cpmvec[lval-2][mval];
           }
            // End scope 2

           //
           //.... Prepare the term in brackets in equation (4.4.11)
           //

           std::complex<double> zbrackets = zpart1 - zpart2;

           //

           std::complex<double> zfactor;

           {
            int const itemp = lval - mval;

            double const dtemp3 = static_cast<double>(itemp);  

            double const dtemp4 = 1.0e+00/dtemp3;

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
   // End of clpmn()

/**
 *  Function: complex_unnormalized_assoc_irregular_legendre()
 *
 *         
 */
void complex_unnormalized_assoc_irregular_legendre(
           unsigned int const mmax, 
           unsigned int const lmax,
           std::complex<double> z,
           std::vector<std::vector<std::complex<double> > > &zcqmvec)
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
      printf("        z = (%10.6f,%10.6f) \n", z.real(), z.imag());
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
  
   std::complex<double> ZMINUS1;

   ZMINUS1.real(-1.0e+00);
   ZMINUS1.imag(0.0e+00);

   std::complex<double> ZERO;  std::complex<double> ZONE;  std::complex<double> ZTWO;

   ZERO.real(0.0e+00);         ZONE.real(1.0e+00);         ZTWO.real(2.0e+00);
   ZERO.imag(0.0e+00);         ZONE.imag(0.0e+00);         ZTWO.imag(0.0e+00);

   //

   std::complex<double> ZHALF;

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

   double const drealabs = z.real();
   double const dimagabs = z.imag();

   bool const z_is_real_and_abs_is_1 = (drealabs == 1.0D+00) && 
                                       (dimagabs == 0.0D+00);

   //

   double const zabs = std::abs(z);

   bool const z_inside_unit_circle  = zabs < 1.0e+00;

   bool const z_outside_unit_circle = zabs >= 1.0e+00;

   //

   if(zdebug)
     {
      printf("      Absolute value of z = %13.6f ", zabs);
      printf("\n\n");

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

   std::complex<double> zls_factor; 

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

   std::complex<double> const zsqd = z * z;

   std::complex<double> const zone_minus_zsqd = ZONE - zsqd;

   std::complex<double> const zs = zls_factor * zone_minus_zsqd;

   std::complex<double> const zq = std::sqrt(zs);  

   //
   //---- Computation of zcq0 
   //
   //      (1/2) * log (s * (z + 1) / (z - 1))
   //
   //     The factor s represent zls_factor above.
   //

   std::complex<double> zcq0; 

   {
    std::complex<double> z_plus_1  = z + ZONE; 

    std::complex<double> z_minus_1 = z - ZONE;

    std::complex<double> ztemp1 = z_plus_1 / z_minus_1;

    std::complex<double> ztemp2 = zls_factor * ztemp1;

    std::complex<double> ztemp3 = std::log(ztemp2); 

    zcq0 = ZHALF * ztemp3;
   }
    // End of scope for computation of zcq0

   if(zdebug)
     {
      printf("      zcq0 = (%13.6e,%13.6e) ", zcq0.real(), zcq0.imag()); 
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

      std::complex<double> zcqf2 = ZERO;
      std::complex<double> zcqf1 = ZONE;
      std::complex<double> zcqf0 = ZERO;

      for(int k = km; k>= 0; k--)
         {
          //
          //.... Ok, build the left hand side term within the bracket
          //

          int const itwok = 2 * k;

          int const itwok_plus_three = itwok + 3;

          double const dtwok_plus_three = static_cast<double>(itwok_plus_three);

          std::complex<double> ztwok_plus_three;

          ztwok_plus_three.real(dtwok_plus_three); 
          ztwok_plus_three.imag(0.0e+00);

          std::complex<double> ztemp1 = ztwok_plus_three * z;

          std::complex<double> znumerator1 = ztemp1 * zcqf1; 

          //
          //.... Ok, build the right hand side term within the bracket
          //

          int const k_plus_two     = k + 2;

          double const dk_plus_two = static_cast<double>(k_plus_two);

          std::complex<double> zkplus_two;

          zkplus_two.real(dk_plus_two);
          zkplus_two.imag(0.0e+00);

          std::complex<double> znumerator2 = zkplus_two * zcqf2; 

          //
          //.... Compute the bracket in equation (4.5.10)
          //
          //     Its the numerator of the fraction
          //

          std::complex<double> znumerator = znumerator1 - znumerator2;

          //
          //.... Compute the denominator of equation (4.5.10)
          //

          int const k_plus_one = k + 1;

          double const dk_plus_one = static_cast<double>(k_plus_one);

          std::complex<double> zk_plus_one;

          zk_plus_one.real(dk_plus_one);
          zk_plus_one.imag(0.0e+00);

          std::complex<double> zdenominator = zk_plus_one;

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

      std::complex<double> zfactor = zcq0 / zcqf0; 

      for(int k=0; k<=lmax; ++k)
         {
          std::complex<double> ztemp = zcqmvec[k][0] * zfactor; 

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

          double const dtwok_plus_three = static_cast<double>(itwok_plus_three);

          std::complex<double> ztwok_plus_three;

          ztwok_plus_three.real(dtwok_plus_three); 
          ztwok_plus_three.imag(0.0e+00);

          std::complex<double> ztemp1 = ztwok_plus_three * z;

          std::complex<double> znumerator1 = ztemp1 * zcqf1; 

          //
          //.... Ok, build the righthand side term within the bracket
          //


          int const k_plus_one     = k + 1;

          double const dk_plus_one = static_cast<double>(k_plus_one);

          std::complex<double> zkplus_one;

          zkplus_one.real(dk_plus_one);
          zkplus_one.imag(0.0e+00);

          std::complex<double> znumerator2 = zkplus_one * zcqf2; 

          //
          //.... Compute the bracket in equation (4.5.10)
          //
          //     Its the numerator of the fraction
          //

          std::complex<double> znumerator = znumerator1 - znumerator2;

          //
          //.... Compute the denominator of equation (4.5.10)
          //

          int const k_plus_two = k + 2;

          double const dk_plus_two = static_cast<double>(k_plus_two);

          std::complex<double> zk_plus_two;

          zk_plus_two.real(dk_plus_two);
          zk_plus_two.imag(0.0e+00);

          std::complex<double> zdenominator = zk_plus_two;

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

      std::complex<double> zfactor1_temp = ZMINUS1 / zq;
 
      std::complex<double> zfactor1 = zfactor1_temp / zcqf0; 

      for(int k=0; k<=lmax; ++k)
         {
          std::complex<double> ztemp = zcqmvec[k][1] * zfactor1; 

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
          std::complex<double> zcq0 = zcqmvec[lval][0];
          std::complex<double> zcq1 = zcqmvec[lval][1];

          for(int mval=2; mval<=mmax; ++mval)
             {
              //
              //.... First term in the recursion of (4.4.10)
              //

              std::complex<double> zterm_first = ZERO;

              {
               double const dmval_plus_1 = static_cast<double>(mval - 1);

               double const dtemp11      = -2.0e+00 * dmval_plus_1;  

               std::complex<double> ztemp11;

               ztemp11.real(dtemp11); ztemp11.imag(0.0e+00);

               std::complex<double> ztemp12 = ztemp11 * z;

               std::complex<double> ztemp13 = ztemp12 / zq;

               zterm_first  = ztemp13 * zcq1;
              }
               // End of scoping for computation of first term in (4.4.10)

              //
              //.... Second term in the recursion of (4.4.10)
              //

              std::complex<double> zterm_second = ZERO;

              {
               double const dlval_plus_mval_minus_1 
                              = static_cast<double>(lval + mval + 1);

               double const dlval_minus_mval_plus_2
                              = static_cast<double>(lval - mval + 2);

               double const dtemp21 = dlval_plus_mval_minus_1 * 
                                        dlval_minus_mval_plus_2;

               std::complex<double> ztemp21;

               ztemp21.real(dtemp21); ztemp21.imag(0.0e+00);

               zterm_second  = ztemp21 * zcq0;
              }
               // End of scoping for computation of second term in (4.4.10)

              //
              //.... Compute the recursion and store 
              //
              //     Note that its a subtraction
              //

              std::complex<double> zcqf = zterm_first - zterm_second;

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

//************************************************************************
//************************************************************************
//
//   End of file 
//
//************************************************************************
//************************************************************************
