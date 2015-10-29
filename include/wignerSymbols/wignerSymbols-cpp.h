/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the LGPL was not -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/lgpl.html.              -/
 ********************************************************/

#ifndef WIGNER_SYMBOLS_CPP_H
#define WIGNER_SYMBOLS_CPP_H

/** \file wignerSymbols-cpp.h
 *
 * 	\author Joey Dumont <joey.dumont@gmail.com>
 *
 * 	\since 2013-08-16
 *
 * 	\brief Defines utility functions for the evaluation of Wigner-3j and -6j symbols.
 *
 * We compute the Wigner-3j and -6j symbols
 * f(L1) = (  L1   L2 L3)
 *         (-M2-M3 M2 M3)
 * for all allowed values of L1, the other parameters
 * being held fixed. The algorithm is based on the work
 * by Schulten and Gordon.
 * 	K. Schulten, "Exact recursive evaluation of 3j- and 6j-coefficients for quantum-mechanical coupling of angular momenta," 
 *		J. Math. Phys. 16, 1961 (1975).
 *  K. Schulten and R. G. Gordon, "Recursive evaluation of 3j and 6j coefficients,"
 *		Comput. Phys. Commun. 11, 269–278 (1976).
 */

#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <iostream>

namespace WignerSymbols {

/*! @name Evaluation of Wigner-3j and -6j symbols.
 * We implement Schulten's algorithm in C++.
 */
///@{
std::vector<double> wigner3j(double l2, double l3,
						double m1, double m2, double m3);

double wigner3j(double l1, double l2, double l3,
					double m1, double m2, double m3);

double wigner3j_auxA(double l1, double l2, double l3,
						double m1, double m2, double m3);
double wigner3j_auxB(double l1, double l2, double l3,
						double m1, double m2, double m3);

std::vector<double> wigner6j(double l2, double l3,
						double l4, double l5, double l6);

double wigner6j(double l1, double l2, double l3,
					double l4, double l5, double l6);

double wigner6j_auxA(double l1, double l2, double l3,
						double l4, double l5, double l6);

double wigner6j_auxB(double l1, double l2, double l3,
						double l4, double l5, double l6);

/*! Computes the Clebsch-Gordan coefficient by relating it to the
 * Wigner 3j symbol. It sometimes eases the notation to use the
 * Clebsch-Gordan coefficients directly. */
inline double clebschGordan(double l1, double l2, double l3,
							double m1, double m2, double m3)
{
	// We simply compute it via the 3j symbol.
	return (pow(-1.0,l1-l2+m3)*sqrt(2.0*l3+1.0)*wigner3j(l1,l2,l3,m1,m2,-m3));
}
}

#endif  // WIGNER_SYMBOLS_CPP_H
