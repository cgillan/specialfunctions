//************************************************************************
//************************************************************************
//
//   associated_Legendre_functions_real.hxx  
//
//   This file implesments templates to compute regular and irregular
//   associated Legendre functions for values on the REAL axis
//   excluding:
//
//      i) The branch cut -infinity to -1 inclusive
//
//     ii) The other pole at +1.
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
//   Copyright (c) 2018,2019 Charles J Gillan  
//   All rights reserved
//
//************************************************************************
//************************************************************************

#ifndef _ASSOC_LEGENDRE_TEMPLATED_REAL_INCLUDE_176492_H_

#define _ASSOC_LEGENDRE_TEMPLATED_REAL_INCLUDE_176492_H_  1

#include <stdio.h>
#include <stdlib.h>

#include <typeinfo> 

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits> 

#include <string>

#include <vector>
#include <algorithm>

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
   unnormalised_associated_irregular_Legendre_big_arg(int const mmax, 
                                                      int const lmax,
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

   T const xsquared = x*x;

   T const bracket  = 1.0e+00 - xsquared;

   T const s_fac = std::signbit(bracket) ? -1.0e+00 : 1.0e+00;

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

   T const smallest = 1.0e+04 * std::numeric_limits<T>::min(); 

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

   T CFL = 0.0e+00;

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
       T ftemp = (n2*x*qlm_array[ll+1][0]) + (n1*qlm_array[ll+2][0]);

       qlm_array[ll][0] = ftemp/n3;

       n1 = n1 + 1;
       n2 = n2 - 2;
       n3 = n3 - 1;
      }

   //
   //---- Renormalize using known value of Q_00 
   //
     
   T const log_factor      = std::log(x + 1.0e+00) - std::log( std::abs(x - 1.0e+00) );

   T const half_log_factor = 0.5e+00 * log_factor;

   T const scale_factor    = half_log_factor/qlm_array[0][0];

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
       T const ftemp = (n2*x*qlm_array[ll+1][1]) + (n1*qlm_array[ll+2][1]);

       qlm_array[ll][1] = ftemp / n3;
             
       n1 = n1 + 1; 
       n2 = n2 - 2; 
       n3 = n3 - 1; 
      }

   //
   //---- Renormalize using known value of first member.
   //
          
   T const root     = std::sqrt ( s_fac * bracket );

   T const Q_01     = -1.0e+00/root;

   T const scale_factor_m_1 = Q_01 / qlm_array[0][1];
      
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

   T const factor = std::sqrt( s_fac * bracket );

   T const factor_inv_bracket  = 1.0e+00/factor;

   for(int ll=0; ll<=lmax; ++ll)
      {   
       for(int m = 1; m<mmax; ++m)
          {
           //
           //.... Factor in front of Qlm = -2 m / std::sqrt( sign (1 - x^2 ) )
           //

           T const factor_Qlm    = -2.0e+00 * (static_cast<double>(m)) * factor_inv_bracket;

           T const term1 = factor_Qlm * x * qlm_array[ll][m];

           //
           //.... Factor in front of Ql{m-1} = - ( sign (1 - x^2 ) ) (l + m)(l - m + 1)
           //

           T const factor_Q_lm_1 = -s_fac * ( static_cast<double>( (ll + m) * (ll - m + 1) ) );

           T const term2 = factor_Q_lm_1 * qlm_array[ll][m-1];

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
 *  Function: unnormalized_assoc_irregular_Legendre_small_arg()
 *
 *  Computes the function for |x| < 1
 *
 */
template <typename T> 
   typename std::enable_if<std::is_floating_point<T>::value,void>::type  
     unnormalized_assoc_irregular_Legendre_small_arg(
           unsigned int const mmax, 
           unsigned int const lmax,
           T xarg,
           std::vector<std::vector<T> > &cqmvec)
  {
   bool const zdebug = true;

   std::string const method_name_str = "unnormalized_assoc_irregular_Legendre_small_arg()";

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
         cformat_str = "        xarg = %10.6f  \n";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "        xarg = %10.6Lf \n";
        }

      printf(cformat_str.c_str(), xarg);

      printf("\n\n");
      printf("      **** End of input data");
      printf("\n\n");
     }

   //
   //---- Establish abs(x) as a constant for rest of routine
   //

   T const xabs = std::abs(xarg);

   //
   //---- Let's establish the constants -2.0, -1.0, 0.0, 0.5, 1.0, 2.0 
   //
  
   T const DMINUS2 = -2.0e+00;
   T const DMINUS1 = -1.0e+00;
   T const DZERO   =  0.0e+00;
   T const DHALF   =  0.5e+00;
   T const DONE    =  1.0e+00;
   T const DTWO    =  2.0e+00;

   //
   //---- Initialise the output vectors 
   //

   for(int lval=0; lval<=lmax; ++lval)
      {
       for(int mval = 0; mval<=mmax; ++mval)
          {
           cqmvec[lval][mval] = DZERO;
          } 
      }

   //
   //--- Are we at the point +1 or -1 on the real axis.
   //
   //    These are the poles of the function 
   //
   //    We stop the computation at this point.
   //
   //    The threshold is somewhat arbitary here.
   //

   T const MULTIPLIER = 20.0e+00;

   T const EPSILON    = std::numeric_limits<T>::epsilon();

   T const THRESHOLD  = MULTIPLIER * EPSILON;

   T const deltaabs   = std::abs(1.0e+00 - xabs);
   
   bool const x_at_pole = (deltaabs < THRESHOLD);  

   if( x_at_pole )
     {
      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float)  ) )
        {
         cformat_str = "    Argument  xarg = %10.6f \n";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "    Argument  xarg = %10.6Lf \n";
        }

      printf(cformat_str.c_str(), xarg);

      printf("\n\n");
      printf("     is too close to one of the poles x=-1 or x=+1 on the \n");
      printf("     real axis.  No computation will be performed.");
      printf("\n\n");

      exit(0);
     }

   //
   //---- Ok, "x" is not at the poles so now decide if "z" lies
   //     inside (even on), or outside the unit circle.
   //

   T radius_unit_circle = (1.0e+00 + EPSILON); 

   //

   bool const x_outside_unit_circle = xabs > radius_unit_circle;

   //

   if(zdebug)
     {
      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float) ) )
        {
         cformat_str = "      Absolute value of xarg = %13.6f - unit radius = %23.13e  ";
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "      Absolute value of xarg = %13.6Lf - unit radius = %23.18Le ";
        }

      printf(cformat_str.c_str(), xabs, radius_unit_circle);
      printf("\n\n");

      //

      if(x_outside_unit_circle)
        {
         printf("      Argument, xarg, lies outside the unit circle");
         printf("\n\n");
         printf("      This routine does not work for that region ");
         printf("\n\n");

         return;
        }
      else
        {
         printf("      Argument, xarg, lies on or inside the unit circle");
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

   T const xls_factor = 1.0e+00; 

   //
   //---- Next compute some auxiliary variables 
   //

   T const xsqd = xarg * xarg;

   T const xone_minus_xsqd = DONE - xsqd;

   T const xs = xls_factor * xone_minus_xsqd;

   T const xq = std::sqrt(xs);  

   //
   //---- Computation of zcq0 
   //
   //      (1/2) * log ( abs( (x + 1) / (x - 1) ) )
   //
   //     The factor s represent zls_factor above.
   //

   T const x_plus_one  = xarg + DONE; 

   T const x_minus_one = xarg - DONE;

   T const xtempor1 = x_plus_one / x_minus_one;

   T const xtempor2 = std::abs( xtempor1 );

   T const xtempor3 = std::log(xtempor2); 

   T const xq0      = DHALF * xtempor3;
  
   // 

   if(zdebug)
     {
      T testprint = 0.0e+00;

      std::string cformat_str = " ";

      if(  ( typeid(T) == typeid(double) ) ||
           ( typeid(T) == typeid(float) ) )
        {
         cformat_str = "      xq0 = %13.6e "; 
        }
      else if( typeid(T) == typeid(long double) )
        {
         cformat_str = "      xq0 = %13.6Le "; 
        }

      printf(cformat_str.c_str(), xq0); 
      printf("\n\n");
     }

   //===============================================================
   //
   //   I N S I D E / O N  U N I T  C I R C L E  -  A B S ( Z ) <= 1 
   //
   //===============================================================
   //
   //---- We can apply upwards recursion 
   //
   //     So we start by setting the 
   //
   //        l=0, m=0   
   //        l=1, m=0
   //        l=0, m=1
   //        l=1, m=1
   //
   //     quadruplet to start the upwards recursion.
   //
   //     Remember that the index in the vector is [l][m]
   //

   //
   //----  Q (l=0,m=0)
   //

   cqmvec[0][0] = xq0;

   //
   //----  Q (l=1, m=0)
   //

   {
    T const temp10 = xarg * xq0;

    cqmvec[1][0] = temp10 - DONE;
   }

   //
   //----  Q (l=0, m=1)
   //
      
   cqmvec[0][1] = DMINUS1 / xq;
 
   //
   //---- Q (l=1, m=1)
   //

   {
    T const x_over_one_minus_xsqd = xarg / xone_minus_xsqd;

    T const xbracket = xq0 + x_over_one_minus_xsqd;

    T const xproduct = xq * xbracket;
 
    cqmvec[1][1] = DMINUS1 * xproduct;
   }

   //
   //---- Ok, we apply recursion now to generate others.
   //
   //     We set m=0 and m=1 and work from l=2,..., lmax
   //

   for(int mval=0; mval<=1; ++mval)
      {
       for(int lval=2; lval<=lmax; ++lval)
          {
           //
           //.... Build up the nominator as a floating point number 
           //
           //     Start with the factors depending on l,m in brackets
           //

           int const ibracket_first  = 2*lval - 1;

           int const ibracket_second = lval + mval - 1;

           T const dbracket_first    = static_cast<T>(ibracket_first);

           T const dbracket_second   = static_cast<T>(ibracket_second);

           //

           T const xfirst_term  = dbracket_first * xarg * cqmvec[lval-1][mval];

           T const xsecond_term = dbracket_second * cqmvec[lval-2][mval]; 

           T const xnumerator   = xfirst_term - xsecond_term;

           //
           //.... Build up the denominator as a floating point number
           //

           int const itemp = lval - mval;

           T   const xdenominator = static_cast<T>(itemp);

           //

           cqmvec[lval][mval] = xnumerator / xdenominator;
          }
           // End loop over lval 
      }
       // End loop over mval 

   //
   //---- Complete the recursion by working on 
   //
   //     l = 0 ,..., lmax and m = 2, ..., mmax
   //

   for(int lval=0; lval<=lmax; ++lval)
      {
       for(int mval=2; mval<=mmax; ++mval)
          {
           //
           //.... Build up the brackets complex numbers
           //

           int const ibracketa = mval - 1;

           int const ibracketb = lval + mval - 1;

           int const ibracketc = lval - mval + 2;

           T const dbracketa   = static_cast<T>(ibracketa);

           T const dbracketb   = static_cast<T>(ibracketb);

           T const dbracketc   = static_cast<T>(ibracketc);

           //
           //.... First term 
           //

           T const xfirst_tempa = DMINUS2 * dbracketa;

           T const xfirst_tempb = xfirst_tempa * xarg;

           T const xfirst_tempc = xfirst_tempb * cqmvec[lval][mval-1]; 

           T const xfirst_term  = xfirst_tempc / xq;

           //
           //.... Second term
           //

           T const xsecond_tempa = dbracketb * dbracketc;

           T const xsecond_tempb = xsecond_tempa * cqmvec[lval][mval-2];

           T const xsecond_term  = xsecond_tempb * xls_factor;

           //

           cqmvec[lval][mval] = xfirst_term - xsecond_term;
          }
           // End loop over mval 
      }
       // End loop over lval

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
   // End of generation of the Q(l,m) values() |x| < 1

/**
 *  Function: unnormalized_assoc_irregular_Legendre()
 *
 *  This function is a gateway which switches to the appropriate function
 *  for the range of x.
 *
 *      i)  |x| < 1
 *     ii)   x  > 1
 *
 *  This is because different computational methods are needed in 
 *  each domain
 *
 */
template <typename T> 
   typename
      std::enable_if< std::is_floating_point<T>::value, 
                      void >::type
   unnormalised_associated_irregular_Legendre(int const mmax, 
                                              int const lmax,
                                              T   const x,
                                              std::vector<std::vector<T> > &qlm_array)
  {
   //
   //---- Initialise the output vectors 
   //

   for(int lval=0; lval<=lmax; ++lval)
      {
       for(int mval = 0; mval<=mmax; ++mval)
          {
           qlm_array[lval][mval] = 0.0e+00;
          } 
      }

   //
   //---- We divide the x-axis including the branch cut and the 
   //     pole at x = 1 into intervals as follows:
   //
   //     [-infty,-1+THRESHOLD ], Inside, [1-THRESHOLD, 1+THRESHOLD], Outside
   //
   //     We do not calculate for x in the middle interval around +1, nor
   //     on the branch cut
   //

   T const MULTIPLIER = 10.0e+00;

   T const EPSILON    = std::numeric_limits<T>::epsilon();

   T const THRESHOLD  = MULTIPLIER * EPSILON;

   //
   //---- Classify X into one of the intervals and when appropriate
   //     call the calculation routine.
   //
   
   if( x < (-1.0e+00 + THRESHOLD) ) 
     {
      //
      //.... On branch cut 
      //

      ;  
     }
   else if( (x > (-1.0e+00 + THRESHOLD) ) && (x < (1.0e+00 - THRESHOLD) ) )
     {
      //
      //.... x in [-1,+1]
      //

      unnormalized_assoc_irregular_Legendre_small_arg(mmax,lmax,x,qlm_array);
     }
   else if( x > (1.0e00 + THRESHOLD) )
     {
      //
      //.... x > 1 (to right of pole)
      //

      unnormalised_associated_irregular_Legendre_big_arg(mmax,lmax,x,qlm_array);
     }
   else
     {
      //
      //.... x lies near pole at 1
      //

      ;
     }
      // End of selection on value of x

   //

   return;
  }
   // End of unnormalised_associated_irregular_Legendre()

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

#endif // For #ifdef _ASSOC_LEGENDRE_TEMPLATED_REAL_INCLUDE_176492_H_

//************************************************************************
//************************************************************************
//
//   End of file 
//
//************************************************************************
//************************************************************************
