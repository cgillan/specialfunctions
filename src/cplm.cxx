//*********************************************************
//*********************************************************
//
//   File: cplm.cxx
//
//   Developed from an initial Fortran 77 version by 
//   Bell and McLaughlin (1983). The original version was 
//   not published but was used in the work for their 
//   latetr in J Phys B on H-H collisions.
//
//   NB: Need --std=c++11 (or similar) on 
//       g++ otherwise this fails 
//
//   Copyright (c) 2018 Charles J Gillan  
//   All rights reserved
//
//*********************************************************
//*********************************************************

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <cmath>
#include <complex>

#include <vector>
#include <map>

double factorial(unsigned int const n);

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

//**********************************************************************
//                                                                      
//     Calculates the regular associated legendre function  P  ( z )    
//                                                           l,m       
//     l,m assumed to be integer                            
//     z   complex and |z| > 1.0 
//                     
//     Developed from an initial version by Bell and McLaughlin (1983)
//                                                 
//**********************************************************************
void cplm(unsigned int const l, 
          unsigned int const m,
          std::complex<double> const &z,
          std::complex<double> &zretVal)                                  
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

   double const  DHALF = 0.5e+00;
   double const  DONE  = 1.0e+00;
   double const  DTWO  = 2.0e+00;

   //
   //.... Local boolean variables  
   //

   bool const zdebug = false;

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

   double const denominator = xtwo_power_l * factorial(l) * factorial(lm);                  

   double       wlmro       = factorial(itwo_l) / denominator;                  

   std::complex<double> zwlmro;

   zwlmro.real(wlmro);  zwlmro.imag(0.0e+00);

   //

   std::complex<double> zpow = power_complx_to_unsigned_int<double,false>(z,lm);

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

       std::complex<double> const zpowterm = power_complx_to_unsigned_int<double,false>(z,lmir);

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
   // End function cplm 

//********************************************************************** 
//
//     Compute factorial of "n"
//
//********************************************************************** 
double factorial(unsigned int const n)
  { 
   double res;

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
      res = 1.0e+00;
     }
   else if(n == 1)
     {
      res = 1.0e+00;
     }
   else if(n > 1)
     {
      int const n1 = n - 1;

      res = ( static_cast<double>(n) ) * factorial(n1);
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

      os << "          Result: " << res 
         << "\n\n" 
         << "          <<<<< Completed: factorial(unsigned int)"
         << "\n\n";

      std::cout << os.str() << "\n";
     }

   return  res;
  }
   // End function factorial

/**
 *   main()
 *
 */
int main(int argc, char **argv)
  {
   //
   //---- Prepare vector of arguments 
   //

   std::vector<std::complex<double> > zarg_vec;

   {
    std::complex<double> zarg;
 
    zarg.real(1.5e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(2.0e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(2.5e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(3.0e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(3.5e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(4.0e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(4.5e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(5.0e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(5.5e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(6.0e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(6.5e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);

    zarg.real(7.0e+00); zarg.imag(0.0e+00);   zarg_vec.push_back(zarg);
   }

   //
   //---- Following are selected values (l,m) from the tables of Zhang 
   //     and Jin (1996)
   //
   //     See pages 118 and following of their book.
   //

   std::vector<int> l_vec_test_p { 1, 2, 3, 10, 2, 3, 4, 10,  3, 4, 5, 10,  4, 5, 6, 10 };
   std::vector<int> m_vec_test_p { 1, 1, 1,  1, 2, 2, 2,  2,  3, 3, 3,  3,  4, 4, 4,  4 };

   //================================================================
   //
   //    L O O P  O V E R   A R G U M E N T S
   //
   //================================================================

   for(int indx=0; indx<zarg_vec.size(); ++indx)
      {
       std::complex<double> const zarg = zarg_vec[indx];

       printf("\n\n     "); 
         for(int icol=2; icol<72;++icol) printf("-"); 

       printf("\n\n");
       printf("     Computed associated Legendre functions of the first kind (regular)");
       printf("\n\n");
       printf("      l    m        Argument (z)            Associated Legendre function    \n");
       printf("     ---  ---  -----------------------   ----------------------------------- \n");

       for(int i=0; i<l_vec_test_p.size(); ++i)
          {
           int const l = l_vec_test_p[i];
           int const m = m_vec_test_p[i];

           if(m > l) continue; // Q_lm can have m > l but not P_lm.

           std::complex<double> zplm;

           cplm(l,m,zarg,zplm);
                                  
           double const xarg = zarg.real();
           double const yarg = zarg.imag();

           double const xplm = zplm.real();
           double const yplm = zplm.imag();

           printf("     %3d  %3d    (%8.4f,%8.4f)     (%16.9e,%16.9e) \n", 
                  l, m, xarg, yarg, xplm, yplm);
          }
      }
       // End of loop over z arguments

   printf("\n\n");

   //
   //---- End of program 
   //

   exit(0);
  }

//********************************************************************** 
//********************************************************************** 
//
//   End of file 
//
//********************************************************************** 
//********************************************************************** 
