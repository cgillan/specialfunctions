//************************************************************************
//************************************************************************
//
//   File: test3.cxx
//
//   Computes P_lm over a square of values in the complex plane,
//   surrounding the origin subject to the limitations:
//
//      (1) avoid the poles at -1 and +1 on the real axis.
//      (2) avoid the branch cut Re(z) \in [\infty,-1]
//
//   Formats the output in columms suitable for plotting using GNUplot
//
//   The user can edit and set the "l" and "m" values of the function
//   that is to be plotted.
//
//--+----1----+----2----+----3----+----4----+----5----+----6----+----7----
//
//   Version history:
//
//      Date      Modifier           Comments
//   ----------   --------   --------------------------------
//   30-Jul-2018   CJG        Initial version
//
//
//   Copyright 2018 (c) Charles J Gillan
//   All rights reserved
//
//************************************************************************
//************************************************************************

/*.................................................*/
/*                                                 */
/*  Includes from standard C++                     */
/*                                                 */
/*.................................................*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <cmath>
#include <chrono>
#include <tuple>

/*.................................................*/
/*                                                 */
/*  Includes for associated legendre polynomials   */
/*                                                 */
/*.................................................*/

#include "associated_legendre_functions.hxx"

#ifdef SKIP_PRINTING_RESULTS
const bool G_PRINT_RESULTS = false;
#else
const bool G_PRINT_RESULTS = true;
#endif

/**
 *  Main program
 *
 */
int main(int argc, char **argv)
  {
   double const xdelta = 0.01e+00;
   double const ydelta = 0.01e+00;

   unsigned int NX = 500;
   unsigned int NY = 500;

   // read NX / NY as parameters
   if (argc >= 3) {
      NX = std::atoi(argv[1]);
      NY = std::atoi(argv[2]);
   }

   if (NX == 0) NX = 500;
   if (NY == 0) NY = 500;

   std::cout << "Running benchmark for NX=" << NX << " NY=" << NY << "\n";

   //
   //---- We compute triangle in "l" and "m" defined by
   //
   //        l=[0,...,LMAX] and m=[0,l] FORALL "l"
   //

   unsigned int const LMAX = 5;

   //

   std::ostringstream os;

   //

   typedef struct plot_point_struct_t
              {
               unsigned int l;              // Holds "order "
               unsigned int m;              //  ""   "degree"
               std::complex<double> zarg;   //  ""   zarg at which to compute function
               std::complex<double> zfunc;  //  ""   vaue of function when computed
              }
               plot_point_struct_t;

   //

   std::vector<plot_point_struct_t> points_vec;
   points_vec.reserve((2*NX) * (2*NY));

   //===========================================================
   //
   //   D E F I N E   G R I D   P O I N T S
   //
   //===========================================================
   //
   //---- Run around the square and define points to be computed
   //

   int lower_x_bound = -1 * NX;
   int upper_x_bound = NX;

   for(int ix=lower_x_bound; ix<upper_x_bound; ++ix)
      {
       double const xreal = static_cast<double>(ix) * xdelta;

       int lower_y_bound = -1 * NY;
       int upper_y_bound = NY;

       for(int iy=lower_y_bound; iy<upper_y_bound; ++iy)
          {
           double const yimag = static_cast<double>(iy) * ydelta;

           //
           //---- Test the restrictions
           //
           //     Branch cut is on [-\infty,-1]
           //
           //        By definition this also avoids the pole at z = -1.
           //
           //     Pole on postive real axis at z = 1
           //

           bool const b_on_branch_cut = (xreal           < -0.99e+00) &&
                                        (std::abs(yimag) < 0.01e+00);

           bool const b_near_positive_pole = (std::abs(xreal - 1.0e+00) < 0.01e+00 ) &&
                                             (std::abs(yimag)           < 0.01e+00);

           //
           //---- Apply the restrictions
           //

           if(b_on_branch_cut || b_near_positive_pole)
             {
              ; // Do nothing - cleaner code than a "break"
             }
           else
             {
              //
              //.... Ok, we build a struct recording this point
              //     and store it
              //
              //     Do this fo all "l" and "m" required
              //

              for(unsigned int l=0; l<=LMAX; ++l)
                 {
                  for(unsigned int m=0; m<=l; ++m)
                     {
                      plot_point_struct_t  point;

                      point.l = l;
                      point.m = m;

                      point.zarg.real(xreal);
                      point.zarg.imag(yimag);

                      point.zfunc.real(0.0e+00);
                      point.zfunc.imag(0.0e+00);

                      //

                      points_vec.push_back(point);
                     }
                      // End loop over "m" values for each "l"
                 }
                  // End loop on l values up to LMAX
             }
              // End of if statement applying branch/pole restrictions
          }
           // End of for loop on "iy" on imaginary axis
      }
       // End of for loop on "ix" on real axis

   //
   //---- Tell the job log hoiw many point we are going to compute
   //

   os.str(""); os.clear();

   os << "\n\n"
      << "   Highest angular momentum - LMAX = " << LMAX
      << "\n\n"
      << "   Number of points at which to compute associated Legendre function = "
      << std::setw(5) << points_vec.size();

   std::cout << os.str() << "\n\n";

   os.str(""); os.clear();

   //===========================================================
   //
   //   C O M P U T E   F U N C T I O N
   //
   //===========================================================
   //
   //---- For all required grid points compute the function
   //
   //     Each point is independent, so possible parallelism.
   //

   auto const start = std::chrono::high_resolution_clock::now();

   //

  // for(auto &point_ref : points_vec)
   #pragma omp parallel for
   for(std::size_t i = 0; i<points_vec.size(); i++) {
      auto& point_ref = points_vec[i];
      {
       cplm(point_ref.l,
            point_ref.m,
            point_ref.zarg,
            point_ref.zfunc);
      }
   }
       // End of loop over points in the tuple vector.

   //

   auto const finish  = std::chrono::high_resolution_clock::now();

   auto const elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count();

   //

   os.str(""); os.clear();

   os << "   Elapsed time for computation of all functions = " << elapsed << " (milliseconds) ";
   os << "\n";

   std::cout << os.str() << "\n";

   if (!G_PRINT_RESULTS) return EXIT_SUCCESS;

   //===========================================================
   //
   //   F O R M A T   F O R   P L O T T I N G
   //
   //===========================================================
   //
   //---- This can be tuned as needed
   //

   os.str(""); os.clear();

   os << "    Index     l   m        Argument                  |arg|      Associated Legendre Func      |func|      \n";
   os << "   --------  --  --  ---------------------------  -----------  ---------------------------  --------------  ";

   std::cout << os.str() << "\n";

   os.str(""); os.clear();

   //

   unsigned int indx = 0;

   for(auto &point_ref : points_vec)
      {
       os.str(""); os.clear();

       std::ios_base::fmtflags old_settings = std::cout.flags();

       std::cout.flags(old_settings);

       //

       os << "  " << std::setw(9) << indx
          << "  " << std::setw(2) << point_ref.l
          << "  " << std::setw(2) << point_ref.m;

       os << "  ("
          << std::right <<  std::scientific << std::setw(11) << std::setprecision(3)
          << point_ref.zarg.real()
          << ", "
          << std::right <<  std::scientific << std::setw(11) << std::setprecision(3)
          << point_ref.zarg.imag()
          << " ) "
          << " ";

       os << std::right << std::scientific << std::setw(11) << std::setprecision(3)
          << std::abs(point_ref.zarg)
          << " | ";

       os << "("
          << std::right <<  std::scientific << std::setw(11) << std::setprecision(3)
          << point_ref.zfunc.real()
          << ", "
          << std::right <<  std::scientific << std::setw(11) << std::setprecision(3)
          << point_ref.zfunc.imag()
          << " ) "
          << "   ";

       os << std::right << std::scientific << std::setw(11) << std::setprecision(3)
          << std::abs(point_ref.zfunc)
          << "   ";

       std::cout << os.str() << "\n";

       std::cout.flags(old_settings);

       //

       ++indx;
      }
       // End of for loop printing points

   //

   os.str(""); os.clear();

   os << "\n\n"
      << "     **** End of computation of associated Legendre functions ";

   std::cout << os.str() << "\n\n";

   os.str(""); os.clear();

   //
   return 0;
  }

//************************************************************************
//************************************************************************
//
//   End of file
//
//************************************************************************
//************************************************************************
