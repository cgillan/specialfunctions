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
   printf("     Reporting status of floating point exception flags.");
   printf("\n\n");

   if(std::fetestexcept(FE_ALL_EXCEPT))
     {
      if(std::fetestexcept(FE_DIVBYZERO)) 
        {
         printf("     **** Division by zero reported");
         printf("\n");
        }

      if(std::fetestexcept(FE_INEXACT)) 
        {
         printf("     **** Inexact computation ");
         printf("\n");
        }

      if(std::fetestexcept(FE_INVALID)) 
        {
         printf("     **** Invalid computation ");
         printf("\n");
        }

      if(std::fetestexcept(FE_OVERFLOW)) 
        {
         printf("     **** Overflow reported ");
         printf("\n");
        }

      if(std::fetestexcept(FE_UNDERFLOW)) 
        {
         printf("     **** Underflow reported ");
         printf("\n");
        }
     }
   else
     {
      printf("     No floating point exception flags are set.");
      printf("\n");
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

