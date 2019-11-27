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
  
#include "monitor_fl_pt_exceptions.hxx"

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
void monitor_fl_pt_exceptions()
  {
   printf("\n\n");

   if(std::fetestexcept(FE_ALL_EXCEPT))
     {
      if(std::fetestexcept(FE_DIVBYZERO)) 
        {
         printf("     **** Floating point division by zero reported");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_INEXACT)) 
        {
         printf("     **** Floating point inexact computation reported");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_INVALID)) 
        {
         printf("     **** Floating point invalid computation reported");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_OVERFLOW)) 
        {
         printf("     **** Floating point overflow reported ");
         printf("\n\n");
        }

      if(std::fetestexcept(FE_UNDERFLOW)) 
        {
         printf("     **** Floating point underflow reported ");
         printf("\n\n");
        }
     }
   else
     {
      printf("     **** NO floating point exception flags were reported.");
      printf("\n\n");
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

