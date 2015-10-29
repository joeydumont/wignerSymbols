/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the LGPL was not -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/lgpl.html.              -/
 ********************************************************/

#ifndef WIGNER_SYMBOLS_FORTRAN_H
#define WIGNER_SYMBOLS_FORTRAN_H

/** \file wignerSymbols-fortran.h
 *
 * \author Joey Dumont <joey.dumont@gmail.com>
 *
 * \since 2013-08-16
 *
 * \brief Defines utility functions for the evaluation of Wigner-3j and -6j symbols.
 *
 * We use modified SLATEC (http://netlib.org/slatec) Fortran subroutines to compute
 * Wigner-3j and -6j symbols. We modified the subroutines so that they do not depend
 * d1mach, r1mach or i1mach as these are obsolete routines. They have been replaced
 * by intrinsics such as huge(), tiny(), epsilon() and spacing(), which are guaranteed
 * to work. The file wignerSymbols.f contain the subroutines and their
 * dependencies.
 *
 * We rely on the ISO C Binding to bind the Fortran subroutines to C++. This method
 * of working insures that proper type casts are performed, among other things. We
 * then use the same method as we did before (extern "C").
 *
 */

#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <iostream>

namespace WignerSymbols {
extern "C"
{
  extern void drc3jj_wrap(double,double,double,double,double*,double*,double*,int,int*);
  extern void drc6j_wrap(double,double,double,double,double,double*,double*,double*,int,int*);
}

/*! Compute a string of Wigner-3j symbols for given l2,l3,m1,m2,m3. */
std::vector<double> wigner3j_f(double l2, double l3, double m1, double m2, double m3);

/*! Computes the Wigner-3j symbol for given l1,l2,l3,m1,m2,m3. We
 * explicitly enforce the selection rules. */
double wigner3j_f(double l1, double l2, double l3, double m1, double m2, double m3);

/*! Computes the Clebsch-Gordan coefficient by relating it to the
 * Wigner 3j symbol. It sometimes eases the notation to use the
 * Clebsch-Gordan coefficients directly. */
inline double clebschGordan_f(double l1, double l2, double l3, double m1, double m2, double m3)
{
  // We simply compute it via the 3j symbol.
  return (pow(-1.0,l1-l2+m3)*sqrt(2.0*l3+1.0)*wigner3j_f(l1,l2,l3,m1,m2,-m3));
}

/*! Computes the Wigner-6j symbol for given, l1, l2, l3, l4, l5, l6.
 * We explicitly enforce the selection rules. */
double wigner6j_f(double l1, double l2, double l3, double l4, double l5, double l6);
}

/*! Computes a string of Wigner-6j symbols for given l2, l3, l4, l5, l6. */
std::vector<double> wigner6j_f(double l2, double l3, double l4, double l5, double l6);

#endif // WIGNER_SYMBOLS_FORTRAN_H
