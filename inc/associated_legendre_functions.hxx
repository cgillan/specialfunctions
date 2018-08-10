//************************************************************************
//************************************************************************
//
//   associated_legendre_templates.hxx  
//
//   Notes: 
//
//     The template structure expression
//
//         std::is_floating_point<T>::value 
//
//     equals "true" if T is float,double,long double 
//     otherwise false.
//
//   std::remove_const<const int> --> int
//
//   NB: Need --std=c++0x (or similar) on 
//       g++ otherwise this fails 
//
//   Copyright (c) 2018 Charles J Gillan  
//   All rights reserved
//
//************************************************************************
//************************************************************************

#ifndef _ASSOC_LEGENDRE_TEMPLATED_REAL_INCLUDE_176492_H_

#define _ASSOC_LEGENDRE_TEMPLATED_REAL_INCLUDE_176492_H_  1

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

/**
 *  Procedure: real_factorial()
 *
 *  Compute factorial of integer "n"
 *
 *  Return floating point type only.
 *
 */ 
template <typename T> 
    typename std::enable_if<std::is_floating_point<T>::value,void>::type  
       real_factorial(unsigned int const n, T &xfact)
  { 
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

   if(0 == n)
     {
      xfact = 1.0e+00;
     }
   else if(n == 1)
     {
      xfact = 1.0e+00;
     }
   else if(n > 1)
     {
      unsigned int const n1 = n - 1;

      T xfact_1; 

      real_factorial(n1,xfact_1);

      T const xtemp = static_cast<T>(n);

      xfact =  xtemp * xfact_1;
     }
   else
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "          **** Error: factorial(unsigned int)"
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

      os << "          Result: " << xfact 
         << "\n\n" 
         << "          <<<<< Completed: factorial(unsigned int)"
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   return;
  }
   // End function factorial

/**
 *   Procedure: factorial2()
 *
 *   Computes n!! - the double factorial
 *
 *   For an even integer n , the double factorial is 
 *   the product of all even integers less than or 
 *   equal to n but greater than or equal to 2. 
 *
 *   For an odd integer p , the double factorial is 
 *   the product of all odd integers less than or equal 
 *   to p and greater than or equal to 1. 
 *
 *   The double factorial values of 0 and -1 are 
 *   defined as equal to 1. Double factorial values 
 *   for integers less than -1 are not defined.
 * 
 */
template <typename T> 
   typename
      std::enable_if<std::is_floating_point<T>::value,T>::type
   factorial2(unsigned int n)
  {
   T fact2 = 1.0e+00;

   switch(n)
     {
      case  0:
      case  1:
              break;

      case  2:
              fact2 = 2.0e+00;
              break;

      case  3:
              fact2 = 3.0e+00;
              break;

      case  4:
              fact2 = 8.0e+00;
              break;

      default:
              for(int i=n; i>1; i=i-2)
                 {
                  fact2 = fact2 * ( static_cast<T>(i) );
                 }
              break;
     }

   return fact2;
  }
   // End of factorial2()

/**
 *   lentz_cfrac()  
 *
 *   Evaluates continued fractions using the method of 
 *   Lentz and Thompson. This is explained in the book
 *   Numerical Recipes in C - see page 165 and following.
 *
 *   This version is adapted from the Fortran 90 code of 
 *   Schneider et al. (2010) in journal CPC.
 *   
 *   It is useful to remember that in this method, the 
 *   continued fraction is evaluated using rational 
 *   approximations, as follows
 *
 *         C_j = A_j / A_{j-1}
 *  
 *         D_j = B_j / B_{j-1}
 *
 *         f_j = f_{j-1} * C_j * D_j  
 *
 */
template <typename T> 
   typename
      std::enable_if< std::is_floating_point<T>::value, 
                      void >::type
   lentz_cfrac(T *f, T x, int nu, int mu)
  {
   bool const zdebug = false;

   std::ostringstream os;

   //
   //---- Banner header 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     ====> lentz_cfrac() <===="; 
      os << "\n\n";
      os << "     Input data: \n";
      os << "       x  (arg)  = " << x  << "\n";
      os << "       nu (lmax) = " << nu << "\n";
      os << "       mu (m)    = " << mu << "\n";
      os << "\n";
      os << "     **** End of input data ";
      os << "\n\n";

      std::cout << os.str();

      os.str(""); os.clear();
     }

   //
   //---- Begin algorithm 
   //

   *f = std::numeric_limits<T>::min();

   T C    = *f;      // Set C_0 = f_0
   T D    = 0.0e+00; // Set D_0 = 0 

   T a    = 1.0e+00; // Starting value for A
   T b    = 0.0e+00; // Starting value for B

   T n_0  = nu;      // n_0 will be the "l" value; we increase 
                     //     by 1 with each iteration

   int count = 0;    // Count number of iterations around loop

   T looptest = 1.0e+00; // looptest will hold (delta_j  - 1)
                         // and we terminate on it when < epsilon

   while (looptest > std::numeric_limits<T>::epsilon() )              
     {
      count = count + 1;

      b = ( ( n_0 + n_0 + 1.0e+00 ) * x ) / ( n_0 + mu );

      //
      //.... Update D:  D_j = b_j + a_j * D_{j-1} 
      //
      //     Apply Thompson correction if D is zero
      //
      //     Note that we actually want the inverse of this 
      //

      D = b + a * D;   // Update D:  D_j = b_j + a_j D_{j-1} 

      if( D == 0.0e+00 ) D = std::numeric_limits<T>::min();

      D = 1.0e+00 / D;

      //
      //.... Update C: C_j = b_j + a_j/C_{j-1}  
      //

      C = b + a / C;

      if( C == 0.0e+00 ) C = std::numeric_limits<T>::min();

      //
      //.... Compute delta_j = C_j * D_j 
      //
      //     Then        f_j = f_{j-1} * delta_j

      T const delta = C * D;

      (*f) = (*f) * delta;

      looptest = std::abs ( delta - 1.0e+00 );

      a = - ( n_0 - mu + 1.0e+00 ) / ( n_0 + mu );

      n_0 = n_0 + 1.0e+00; // Augment the value of "l" 
     }
      // End of while loop 

   //
   //---- Return point 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     Computed value of continued fraction = " << *f;
      os << "\n\n";
      os << "     **** Completed - lentz_cfrac() "; 
      os << "\n\n";

      std::cout << os.str();

      os.str(""); os.clear();
     }

   return;
  }
   // End of lentz_cfrac()

/**
 *    Method: unnormalised_associated_regular_Legendre()
 *
 *    Computes the unnormalized associated egendre functions P_{lm}
 *    for integer l and integer m.
 *
 *    Limits:   l > 0 and  m in [0,l]
 *
 *    Notes:
 *
 *      This works for x in [-1,1] and in [1,\infty] 
 *
 *      Properties are defined in the reference:  
 *
 *                        http://dlmf.nist.gov/14
 *
 *      The operation of the routine can be explained in 
 *      a visual way as follows. We can draw a half triangle
 *      with "m" index along the top and "l" down the left.
 *      The elements are the P_{lm} (x) function, all at some 
 *      fixed value of argument "x".
 *
 *              <---- m index ----->
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
 *     Given a value "m", say m=5, we start at P00 (top left)
 *     use recursion to generate 
 *
 *           P11, P22, P33, P44, P55. 
 *
 *     Thus we work along the diagonal. 
 *
 *     We can use a different recursion formula to move down
 *     the column from P55. That is we may generate the
 *
 *            P65, P75, P85, P95
 *
 *     Terminating at "lmax" 
 * 
 *     We can exercise this procedure in a loop over "m"
 *     and build a full triangle of P_{lm} values up to
 *     some mmax lmax
 *
 */
template <typename T> 
   typename
      std::enable_if< std::is_floating_point<T>::value, 
                      void >::type
      unnormalised_associated_regular_Legendre(int    const lmax, 
                                               int    const m, 
                                               T      const x,
                                               std::vector<T> &plm_vec)
  {
   bool const zdebug = false;

   std::ostringstream os;

   //
   //---- Banner header 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     ====> unnormalised_associated_regular_Legendre() <===="; 
      os << "\n\n";
      os << "     Input data: \n";
      os << "       lmax  = " << lmax  << "\n";
      os << "       m     = " << m     << "\n";
      os << "       x     = " << x     << "\n";
      os << "\n";
      os << "     **** End of input data ";
      os << "\n\n";

      std::cout << os.str();

      os.str(""); os.clear();
     }

   //
   //---- Exclude the singularities
   //
   //     It is useful to remember that -1,+1 and \infty
   //     are singularities of the associated Legendre 
   //     equation
   //

   //
   //---- Zeroise plm_vec and check size.
   //
   //     NB: Had tried a swap operation to release memory
   //         but this core dumps in "free()" 
   //
   
   if(plm_vec.size() < lmax+1)
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     **** Error in: unnormalised_associated_regular_Legendre() "; 
      os << "\n\n";
      os << "     Input data: \n";
      os << "       lmax  = " << lmax  << "\n";
      os << "       m     = " << m     << "\n";
      os << "       x     = " << x     << "\n";
      os << "\n";
      os << "     Size of storage plm_vec should be (lmax+1) doubles =  " << lmax+1;
      os << "\n\n";
      os << "     But is actually too small at " << plm_vec.size() << " doubles "; 

      std::cout << os.str() << "\n";
     }
   else
     {
      for(int ikkk=0; ikkk<plm_vec.size(); ++ikkk)
         {
          plm_vec[ikkk] = 0.0e+00;
         }
     }

   //
   //---- Decide whether we are in [-1,+1] or outside 
   //
   //     Being inside the range is known as being "on the cut"
   //
   //        On  the cut we have s_fac =  1
   //        Off the cut we have s_fac = -1
   //

   T const xsquared = x*x;

   T const bracket  = 1.0e+00 - xsquared;

   T const s_fac = std::signbit(bracket) ? -1.0e+00 : 1.0e+00;

   //
   //---- Compute P_mm 
   //
   //     Recurse upwards to P_{m,m} from P_{0,0} = 1 
   //
   //     P_mm = -s_fac*(2*m - 1)*sqrt(s_fac*(1 - x*x))*P_(m-1)(m-1)
   //
   //     This is explained in the notes above as moving along 
   //     the diagonal.
   //

   T const P_00 = 1.0e+00;

   plm_vec[0] = P_00;

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "     --------------------------------------------- ";
      os << "\n\n";
      os << "     Starting recursion on \"m\" ";
      os << "\n\n";
      os << "       P(l=0,m=0) = " << P_00 << " at x = " << x;
      os << "\n\n";
   
      std::cout << os.str();

      os.str(""); os.clear();
     }
   
   //

   T ftemp = P_00;

   if(m > 0)
     {
      T const minus_s_fac = -s_fac;
 
      T const root        = sqrt(s_fac * bracket);

      for(int mtemp = 1; mtemp<=m; ++mtemp)
         {
          T const factorm = (double) (2*mtemp - 1); 

          ftemp   = minus_s_fac * factorm * root * ftemp;

          if(zdebug)
            {
             os.str(""); os.clear();

             os << "     --------------------------------------------- ";
             os << "\n\n";
             os << "     Recursion upwards on \"m\" ";
             os << "\n\n";
             os << "       P(l=" << mtemp << ",m=" << mtemp << ") = " << ftemp << " at x = " << x;
             os << "\n\n";
   
             std::cout << os.str();

             os.str(""); os.clear();
            }
         }
     }
      // End of test on m>0 

   //
   //---- Ok, we not can store P_{mm} before we start
   //     to iterature "up" in the column, that is 
   //     in "l" values for fixed "m".
   //

   T P_mm = ftemp;

   plm_vec[m]  = ftemp;
 
   if(zdebug)
     {
      os.str(""); os.clear();

      os << "     --------------------------------------------- ";
      os << "\n\n";
      os << "     Recursion over \"m\" is complete from 0 to " << m;
      os << "\n\n";
      os << "     P(l=" << m << ",m=" << m << ") = " << ftemp << " at x = " << x;
      os << "\n\n";
   
      std::cout << os.str();

      os.str(""); os.clear();
     }

   //
   //---- If Lmax equals M then we have no more work to do
   //

   if(lmax == m)
     {
      if(zdebug)
        {
         os.str(""); os.clear();

         os << "     --------------------------------------------- ";
         os << "\n\n";
         os << "     lmax equals m - no more work to do";
         os << "\n\n";
   
         std::cout << os.str();

         os.str(""); os.clear();
        }
     }
   else
   {

   //
   //---- Consider recurrence relation, equation (2) of the paper
   //     Schneider et al (2010). Let "l = m" in that equation and 
   //     consider what it gives us for P_{l,m}(x), then we have:
   //
   //     P_{m+1,m}(x) + (2m + 1)xP_{m,m}(x) - (2m)P_{m-1,m}(x) = 0
   //
   //     but by definition we have P_{m-1,m}=0. This is because we 
   //     have m <= l for each P_{l,m}. We cannot have l=m-1 therefore.
   //
   //     NOTE: The recurrence printed in eqn(2) of the paper has an
   //           error - there is an "x" fractor in the right most
   //           term - it should not be there. Check other sources
   //           for the recurrence relationships.
   //
   //     Rearranging gives for l=m+1
   //
   //       P_{m+1,m}(x) = (2m + 1)xP_{m,m)(x)
   //  
   //     Thus starting from the P_{m,m} that we have calculated above,
   //     then we can compute by recursion P_{l,m} for l = m+1 to
   //     some upper value, Lmax. P_{m+1,m} is a specical case of the 
   //     recursion formula as shown above.
   //
   //     More generally, from eqn (2) of the paper we have 
   //
   //     (m - l - 1) P_{l+1,m}(x) = -(2l + 1)xP_{l,m}(x) 
   //                                       + (l + m)P_{l-1,m}(x)
   //
   //     for l=m+1 etc...
   //
   //---- Next case: l = "m+1" 
   //
   //        P_{m+1,m}(x) = (2m + 1)xP_{m,m)(x)
   //

   T const two_m        = static_cast<T>(2*m);
   T const two_m_plus_1 = static_cast<T>(2*m + 1);

   T P_l_minus_1 = 0.0e+00;
   T P_l         = P_mm;

   T P_next      = two_m_plus_1 * x * P_mm;

   plm_vec[m+1]  = P_next;

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "     --------------------------------------------- ";
      os << "\n\n";
      os << "     First recursion upwards on \"l\" ";
      os << "\n\n";
      os << "       P(l=" << m+1 << ",m=" << m << ") = " << P_next << " at x = " << x;
      os << "\n\n";
   
      std::cout << os.str();

      os.str(""); os.clear();
     }
         
   //
   //---- Full recursion formula from "l=m+2" onwards
   //

   P_l_minus_1 = P_mm;
   P_l         = P_next;

   for(int ltemp=m+2; ltemp<=lmax; ++ltemp)
      {
       int const l = ltemp - 1; // Current "l" value

       T const term1    = ( static_cast<T>(2*l + 1) ) * x * P_l; 

       T const term2    = ( static_cast<T>(l + m)   ) * P_l_minus_1;

       int    const ibracket = m - l - 1;

       T const bracket  = static_cast<T>(ibracket);

       T const P_l_plus_1 = (-term1 + term2) / bracket;

       if(zdebug)
         {
          os.str(""); os.clear();

          os << "     --------------------------------------------- ";
          os << "\n\n";
          os << "     Next recursion upwards on \"l\" ";
          os << "\n\n";
          os << "       P(l=" << ltemp << ",m=" << m << ") = " << P_l_plus_1 << " at x = " << x;
          os << "\n\n";
   
          std::cout << os.str();

          os.str(""); os.clear();
         }

       //
       //.... Store value 
       //

       plm_vec[ltemp] = P_l_plus_1;

       //
       //.... Push down for next iteration 
       //

       P_l_minus_1 = P_l;

       P_l         = P_l_plus_1;
      }
       // End of upwards recursion on "l"

   }

   //
   //---- Return point
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "     Size of plm_vec[] = " << plm_vec.size();
      os << "\n\n";
      os << "     **** Completed - unnormalised_associated_regular_Legendre() "; 
      os << "\n\n";

      std::cout << os.str();

      os.str(""); os.clear();
     }

   return;
  }
   // End of unnormalised_associated_regular_Legendre()

/**
 *            
 *   Computes irregular associated Legendre functions, Q_LM(x) when x > 1.
 *
 *   Uses backward recursion and the modified Miller algorithm.
 *
 *     (L + 1 - M) T_LM = (2*L + 3) z T_(L+1)M - (L - M + 2) T_(L+2)M
 *
 *  Starting at a large value of L set the last value to zero and the
 *  next to last to one.  Then recur downward which is the stable direction.
 *  The T's are proportional to the desired Q functions.  
 *  The proportionality constant is determined by the known value of Q00.
 *  This allows us to compute the Q's for m=0. 
 *
 *  The process is repeated for Q_01
 *
 */
template <typename T> 
   typename
      std::enable_if< std::is_floating_point<T>::value, 
                      void >::type
   unnormalised_associated_irregular_Legendre_big_arg(int const lmax, 
                                                      int const mmax,
                                                      T   const x,
                                                      std::vector<std::vector<T> > &qlm_array)
  {
   bool const zdebug = false; 

   std::ostringstream os;

   //
   //---- Banner header 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     ====> unnormalised_associated_irregular_Legendre_big_arg() <===="; 
      os << "\n\n";
      os << "     Input data: \n";
      os << "       lmax  = " << lmax  << "\n";
      os << "       mmax  = " << mmax  << "\n";
      os << "       x     = " << x     << "\n";
      os << "\n";
      os << "     **** End of input data ";
      os << "\n\n";

      std::cout << os.str();

      os.str(""); os.clear();
     }

   //
   //---- Exclude the singularities
   //
   //     It is useful to remember that -1,+1 and \infty
   //     are singularities of the associated Legendre 
   //     equation
   //

   //
   //---- Decide whether we are in [-1,+1] or outside 
   //
   //     This implementation usign Miller's algorithm is 
   //     only valid if we are off the cut, that is we
   //     are outside of [-1,+1]
   //

   double const xsquared = x*x;

   double const bracket  = 1.0e+00 - xsquared;

   double const s_fac = std::signbit(bracket) ? -1.0e+00 : 1.0e+00;

   if(std::signbit(bracket))
     {
      ;
     }
   else
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     **** Error in: unnormalised_associated_irregular_Legendre_big_arg()"; 
      os << "\n\n";
      os << "     Input x = " << x;
      os << "\n\n";
      os << "     Should lie outside the cut [-1,+1] ";
      os << "\n\n";

      std::cout << os.str();  // Could be std::err ??

      os.str(""); os.clear();

      exit(0);
     }
      // End of test for argument x off the cut [-1,+1]
 
   //

   double const smallest = 1.0e+04 * std::numeric_limits<double>::min(); 

   //
   //---- Dimension qlm_array correctly
   //

   qlm_array.resize(lmax+1);

   for(int i=0; i<qlm_array.size(); ++i)
      {
       qlm_array[i].resize(mmax+1);
      }

   //===========================================================
   //
   //   D O W N W A R D    R E C U R S I O N   F O R   M = 0
   //
   //===========================================================
   //
   //---- Compute continued fraction Q_L_max/Q_(L_max-1) for m=0
   //

   int m = 0;

   double CFL = 0.0e+00;

   lentz_cfrac(&CFL, x, lmax, m);

   qlm_array[lmax][0]   = CFL * smallest;
   qlm_array[lmax-1][0] = smallest;

   //
   //----  Downward recursion for m = 0
   //    
   //       (L + 1 - M) T_LM = (2*L + 3) x T_(L+1)M + (M - L - 2) T_(L+2)M
   //
   //      First element in the downwards iteration is L = (lmax - 2)
   //      so following computes the factors in front of the terms 
   //      for that case. We modify these terms in the loop for the 
   //      "next" iteration of the loop. 
   //
   //      e.g.  (L + 1 - M) => (lmax - 2 + 1 - 0)   = lmax   - 1 (n3)
   //            (2*L + 3)   => (2*(lmax - 2) + 3)   = 2*lmax - 1 (n2)
   //            (M - L - 2) => (0 - (lmax - 2) - 2) = -lmax      (n1)  
   //

   int n1 = m - lmax; 
   int n2 = lmax + lmax - 1;
   int n3 = lmax - 1;

   for(int ll=lmax-2; ll>=0; --ll)
      {
       double ftemp = (n2*x*qlm_array[ll+1][0]) + (n1*qlm_array[ll+2][0]);

       qlm_array[ll][0] = ftemp/n3;

       n1 = n1 + 1;
       n2 = n2 - 2;
       n3 = n3 - 1;
      }

   //
   //---- Renormalize using known value of Q_00 
   //
     
   double const log_factor      = log(x + 1.0e+00) - log( std::abs(x - 1.0e+00) );

   double const half_log_factor = 0.5e+00 * log_factor;

   double const scale_factor    = half_log_factor/qlm_array[0][0];

   //

   for(int ll=0; ll<=lmax; ++ll)
      {
       qlm_array[ll][0] = scale_factor * qlm_array[ll][0];
      }

   //===========================================================
   //
   //   R E P E A T    F O R   M = 1
   //
   //===========================================================

   if(mmax < 1) return;

   //
   //---- Compute continued fraction Q_lmax / Q_(lmax-1) for m=1
   //
     
   m = 1; 

   CFL = 0.0e+00;

   lentz_cfrac(&CFL, x, lmax, m);

   qlm_array[lmax][1]   = CFL * smallest;
   qlm_array[lmax-1][1] = smallest;

   //
   //---- Downward recursion - as above 
   //

   n1 = 1 - lmax;
   n2 = lmax + lmax - 1; 
   n3 = lmax;

   for(int ll=lmax-2; ll>=0; --ll)
      {
       double const ftemp = (n2*x*qlm_array[ll+1][1]) + (n1*qlm_array[ll+2][1]);

       qlm_array[ll][1] = ftemp / n3;
             
       n1 = n1 + 1; 
       n2 = n2 - 2; 
       n3 = n3 - 1; 
      }

   //
   //---- Renormalize using known value of first member.
   //
          
   double const root     = sqrt ( s_fac * bracket );

   double const Q_01     = -1.0e+00/root;

   double const scale_factor_m_1 = Q_01 / qlm_array[0][1];
      
   for(int ll=0; ll<=lmax; ++ll)
      {
       qlm_array[ll][1] = scale_factor_m_1 * qlm_array[ll][1];
      }

   //=========================================================
   //
   //   R E C U R S E   U P W A R D S   I N   M   P E R   L 
   //
   //=========================================================
   //
   //---- For each l value, step up in m
   //

   double const factor = sqrt( s_fac * bracket );

   double const factor_inv_bracket  = 1.0e+00/factor;

   for(int ll=0; ll<=lmax; ++ll)
      {   
       for(int m = 1; m<mmax; ++m)
          {
           //
           //.... Factor in front of Qlm = -2 m / sqrt( sign (1 - x^2 ) )
           //

           double const factor_Qlm    = -2.0e+00 * ((double) m) * factor_inv_bracket;

           double const term1 = factor_Qlm * x * qlm_array[ll][m];

           //
           //.... Factor in front of Ql{m-1} = - ( sign (1 - x^2 ) ) (l + m)(l - m + 1)
           //

           double const factor_Q_lm_1 = -s_fac * ( (double) ( (ll + m) * (ll - m + 1) ) );

           double const term2 = factor_Q_lm_1 * qlm_array[ll][m-1];

           //  
           //.... Apply forward recursion formula
           //

           qlm_array[ll][m+1] = term1 + term2;
          }
           // End of loop over "m"
      }
       // End of loop over "l"

   //
   //---- Return point 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n";
      os << "     **** Completed - unnormalised_associated_irregular_Legendre_big_arg() "; 
      os << "\n\n";

      std::cout << os.str();

      os.str(""); os.clear();
     }

   return;
  }
   // End of unnormalised_associated_irregular_Legendre_big_arg() 

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
template <typename T> 
   std::complex< typename std::enable_if<std::is_floating_point<T>::value,T>::type > 
           power_complx_to_unsigned_int(std::complex<T> const &z, 
                                        unsigned int    const ipow)
  {
   bool const zdebug = false;

   if(zdebug)
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

   if(zdebug)
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
 *                                                                      
 *     Calculates the regular associated legendre function  P  ( z )    
 *                                                           l,m       
 *     l,m are positive integers                            
 *     z   complex and |z| > (1.0 + eps) 
 *                     
 *     This works on the real axis - just set the complex part of the 
 *     number "z" to zero on calling the function.
 *  
 *     Developed from a Fortran 77 version by Bell and McLaughlin (1983)
 *     That version was never published.                                          
 *
 */
template <typename T> 
   typename std::enable_if<std::is_floating_point<T>::value,void>::type 
     cplm(unsigned int const l, 
          unsigned int const m,
          std::complex<T> const z,
          std::complex<T> &zretVal)                                  
  {
   std::complex<T> const ZERO(0.0e+00, 0.0e+00); 
   std::complex<T> const ZONE(1.0e+00, 0.0e+00); 
   std::complex<T> const ZTWO(2.0e+00, 0.0e+00);

   //

   T const  DTHRESHOLD_POLE = 0.01e+00; 

   T const  DHALF = 0.5d+00;
   T const  DONE  = 1.0d+00;
   T const  DTWO  = 2.0d+00;

   //

   bool const zdebug = false;

   //

   std::ostringstream os;

   //
   //---- Debug banner header 
   //

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     >>>> cplm - compute regular associated Legendre functions - complex plane "
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
   //     l should be >= 0   - controlled via arg type - insigned int  
   //     m should be [0,+l] - unsigned int so guarantee to in [0,...]
   //
   //     z should NOT be on the branch cut [-infinity,-1] or near z=1.
   //

   if(m > l)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     ***** Error in cplm - compute associated Legendre functions"
         << "\n\n"
         << "     l is greater than m " 
         << "\n\n"
         << "     l = " << std::setw(3) << l << " m = " << std::setw(3) << m
         << "\n";

      std::cout << os.str() << "\n";

      exit(-1);
     }

   //

   T const realPart = z.real();
   T const imagPart = z.imag();

   bool const b_on_branch_cut = (realPart           < -0.99e+00) && 
                                (std::abs(imagPart) < 0.01e+00);

   bool const b_near_positive_pole = (std::abs(realPart - 1.0e+00) < DTHRESHOLD_POLE ) &&
                                     (std::abs(imagPart)           < 0.01e+00);


   if( b_on_branch_cut || b_near_positive_pole )
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     ***** Error in cplm - compute associated Legendre functions"
         << "\n\n";

      if( b_on_branch_cut )
        {
         os << "     z is on the branch cut [-infinity,-1] "; 
        }
      else if( b_near_positive_pole )
        {
         os << "     z is near the pole at z=1 ";
        }
      else
        {
         os << "     Unknown error ";
        }

      os << "\n\n"
         << "     z = ( "
         << std::scientific << std::setw(13) << std::setprecision(7) 
         << z.real() << " , " << z.imag() << " ) " 
         << "   "
         << "\n";

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
   //---- Fast computation of 2^l 
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

   T const xtwo_power_l = static_cast<T>(itwo_power_l);

   if(zdebug)
     {
      os.str(""); os.clear();

      os << "\n\n"
         << "     l = " << l << ",  2^l = " << itwo_power_l << " and as float = " << xtwo_power_l
         << "\n\n";

      std::cout << os.str() << "\n";

     } 

   //
   //---- Compute the required factorials
   //

   T l_factorial = 0.0e+00;

   real_factorial(l, l_factorial);

   //

   T lm_factorial = 0.0e+00;

   real_factorial(lm, lm_factorial);

   //

   T two_l_factorial = 0.0e+00;

   real_factorial(itwo_l,two_l_factorial); 

   //
   //---- Put the pieces together to compute the first term 
   //

   T const denominator = xtwo_power_l * l_factorial * lm_factorial;   
                
   T       wlmro       = two_l_factorial / denominator;                  

   std::complex<T> zwlmro;

   zwlmro.real(wlmro);  zwlmro.imag(0.0e+00);

   //

   std::complex<T> zpow = power_complx_to_unsigned_int<T>(z,lm);

   std::complex<T> zsum = zwlmro * zpow;

   //std::cout << "     Debug: outside loop - power of z (lm) = " << lm << "\n";

   //
   //---- Now deal with terms ir = 1, 2, ... lm2 - indeed if any 
   //                                  
   //     We compute the following term, "wlmr",  at each value "ir" 
   //
   //            -1 * (l - m - 2*ir + 2) * (l - m -2*ir + 1)
   //           --------------------------------------------- * z^(l - m - 2*ir)
   //                 (2*ir) * ( 2(l - ir) + 1 )
   //
   //     Note the product relationship as we save this to "wlmro" and reuse
   //     at next iteration.
   //
 
   for(int ir=1; ir<=lm2; ++ir)
      {
       int const lmir = lm - 2*ir;                                                 
       int const jlir = l  - ir;                                                  

       T wlmr = -wlmro*(static_cast<T>(lmir) + DTWO)*(static_cast<T>(lmir) + DONE);                       
       
       wlmr = wlmr/(DTWO*static_cast<T>(ir)*(DTWO*static_cast<T>(jlir)+ DONE));                

       //

       std::complex<T> zwlmr;

       zwlmr.real(wlmr); zwlmr.imag(0.0e+00);

       //

       //std::cout << "     Debug: inside  loop - ir = " << ir 
       //          << " power of z (lmir) = " << lmir 
       //          << "\n";

       std::complex<T> const zpowterm = power_complx_to_unsigned_int<T>(z,lmir);

       std::complex<T> const t = zwlmr * zpowterm;                                             

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
      T const xn  = static_cast<T>(m) * DHALF;                                          

      std::complex<T> const zsqd   = z*z;

      std::complex<T> const ztemp  = zsqd - ZONE;

      std::complex<T> const ztemp2 = std::pow(ztemp, xn);

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
   // End function cplm 


/**
 *   Procedure: printRectangularQlmMat()
 *
 *   Prints a matrix of Qlm values stored in rectangular form.
 *
 *   Typically this is a set of Q_lm functions computed at one argument
 *
 *   Rows are 0 to lmax
 *
 *   Cols are 0 to mmax
 * 
 *   Storage assumed std::vector<std::vector<T> > 
 * 
 *   T can be  float, double or long double
 *
 */
template <typename T> void printRectangularQlmMat(int const lmax,
                                                  int const mmax,
                                                  std::vector<std::vector<T> > const &mat)
  {
   std::ostringstream os;

   //
   //---- Header on job log 
   //

   std::ios_base::fmtflags old_settings = std::cout.flags();

   os.str(""); os.clear();

   os << "\n\n"             
      << "     Printing Rectangular matrix of Qlm functions  "
      << "\n\n" 
      << "     Number of rows = " << std::fixed << std::setw(5) << lmax + 1 << "\n" 
      << "     Number of cols = " << std::fixed << std::setw(5) << mmax + 1 
      << "\n\n";
   
   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   std::cout.flags(old_settings);

   //
   //.... Row of "m=x" accross the top line
   //
   //     Each "m" is allocated 6 chars and is set 
   //     in a box size 13 characters.
   //

   old_settings = std::cout.flags();

   os.str(""); os.clear();

   os << "          " ;

   for(int mm=0; mm<=mmax; ++mm)
      {
       os << "      m =" << std::setw(2) << std::fixed << std::right << mm << "     ";
      }

   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   std::cout.flags(old_settings);

   //
   //.... Row of "-----" accross under each "m"
   //

   old_settings = std::cout.flags();

   os.str(""); os.clear();

   os << "          ";

   for(int mm=0; mm<=mmax; ++mm)
     {
      os << " +------------+ ";
     }

   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   std::cout.flags(old_settings);

   //
   //.... For each "l" print out the m values 
   //

   for(int ll=0; ll<=lmax; ++ll)
      {
       os.str(""); os.clear();

       os << "     l=" << std::setw(2) << ll << " ";

       for(int mm=0; mm<=mmax; ++mm)
          {
           os << " " 
              << std::scientific << std::setw(14) << std::setprecision(7) 
              << mat[ll][mm] 
              << " " ;
          }

       std::cout << os.str() << "\n";
      }
  
   std::cout.flags(old_settings);

   //

   os.str(""); os.clear();

   os << "\n\n";
   os << "     **** Completed - printLowerTriangularQlmMat()"; 
   os << "\n\n";

   std::cout << os.str();

   os.str(""); os.clear();

   //
   //--- Return point 
   //

   return;
  }
   // End of printRectangularQlmMat()

/**
 *   Procedure: printLowerTriangularPlmMat()
 *
 *   Prints a matrix of regular Plm values stored in lower triangular form.
 *
 *   Typically this is a set of P_lm functions computed at one argument
 *
 *   Rows from 0 to lmax - hence count = lmax+1
 *
 *   Row 0 has 1 element, row 1 has 2 elements etc...
 * 
 *   Storage assumed std::vector<std::vector<T> > 
 * 
 *   T can be  float, double or long double
 *
 */
template <typename T> void printLowerTriangularPlmMat(std::vector<std::vector<T> > const &triang_mat, 
                                                      int const lmax)
  {
   std::ostringstream os;

   //

   int const mmax = lmax+1;

   //
   //---- Header on job log 
   //

   std::ios_base::fmtflags old_settings = std::cout.flags();

   os.str(""); os.clear();

   os << "\n\n"             
      << "     Printing lower triangular matrix  "
      << "\n\n" 
      << "     Number of rows = " << std::fixed << std::setw(5) << lmax + 1 
      << "\n\n";
   
   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   std::cout.flags(old_settings);

   //
   //.... Row of "m=x" accross the top line
   //
   //     Each "m" is allocated 6 chars and is set 
   //     in a box size 13 characters.
   //

   old_settings = std::cout.flags();

   os.str(""); os.clear();

   //os << "1234567890" ;

   os << "          " ;

   for(int mm=0; mm<mmax; ++mm)
      {
       os << "      m =" << std::setw(2) << std::fixed << std::right << mm << "    ";
      }

   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   std::cout.flags(old_settings);

   //
   //.... Row of "-----" accross under each "m"
   //

   old_settings = std::cout.flags();

   os.str(""); os.clear();

   //os << "1234567890";

   os << "          ";

   for(int mm=0; mm<mmax; ++mm)
     {
      os << " +-----------+ ";
     }

   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   std::cout.flags(old_settings);

   //
   //.... For each "l" print out the m values 
   //

   for(int ll=0; ll<=lmax; ++ll)
      {
       os.str(""); os.clear();

       os << "     l=" << std::setw(2) << ll << " ";

       for(int mm=0; mm<=ll; ++mm)
          {
           os << " " 
              << std::scientific << std::setw(13) << std::setprecision(7) 
              << triang_mat[ll][mm] 
              << " " ;
          }

       std::cout << os.str() << "\n";
      }
  
   std::cout.flags(old_settings);

   //

   os.str(""); os.clear();

   os << "\n\n";
   os << "     **** Completed - printLowerTriangularPlmMat()"; 
   os << "\n\n";

   std::cout << os.str();

   os.str(""); os.clear();

   //
   //--- Return point 
   //

   return;
  }
   // End of printLowerTrianglarPlmMat() 

#endif // For #ifdef _ASSOC_LEGENDRE_TEMPLATED_REAL_INCLUDE_176492_H_

//************************************************************************
//************************************************************************
//
//   End of file 
//
//************************************************************************
//************************************************************************
