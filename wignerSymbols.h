#ifndef WIGNER_SYMBOLS_H
#define WIGNER_SYMBOLS_H

/** \file wignerSymbols.h  
 *
 * 	\author Joey Dumont <joey.dumont@gmail.com>
 *
 * 	\since 2013-08-16
 *
 * 	\brief Defines utility functions for the evaluation of Wigner-3j and -6j symbols. 
 *
 * We use modified SLATEC (http://netlib.org/slatec) Fortran subroutines to compute 
 * Wigner-3j and -6j symbols. We modified the subroutines so that they do not depend 
 * d1mach, r1mach or i1mach as these are obsolete routines. They have been replaced 
 * by intrinsics such as huge(), tiny(), epsilon() and spacing(), which are guaranteed
 * to work. The files wignerSymbols.f90 and utils.f contain the subroutines and their
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
#include <armadillo>

namespace varPhase {
/** @name Fortran Linkage
 * We link the Fortran subroutines to our C++ code. 
 * Because of the ISO C Binding (see file wigner_wrap.f90), 
 * some arguments are passed by value rather than by reference. 
 */
///@{
extern "C"
{
	/*! Wigner-3j symbols. */
	extern void drc3jj_wrap(double,double,double,double,double*,double*,double*,int,int*);

	/*! Wigner-6j symbols. */
	extern void drc6j_wrap(double,double,double,double,double,double*,double*,double*,int,int*);
}
///@}

/*! @name Evaluation of Wigner-3j and -6j symbols. 
 * We implement the Fortran subroutines in C++. 
 */
///@{
arma::colvec c_wigner3j(double l2, double l3,
						double m1, double m2, double m3);

double c_wigner3j(double l1, double l2, double l3,
					double m1, double m2, double m3);

double c_wigner3j_auxA(double l1, double l2, double l3,
						double m1, double m2, double m3);
double c_wigner3j_auxB(double l1, double l2, double l3,
						double m1, double m2, double m3);

arma::colvec c_wigner6j(double l2, double l3,
						double l4, double l5, double l6);

double c_wigner6j(double l1, double l2, double l3,
					double l4, double l5, double l6);

double c_wigner6j_auxA(double l1, double l2, double l3,
						double l4, double l5, double l6);

double c_wigner6j_auxB(double l1, double l2, double l3,
						double l4, double l5, double l6);

/*! Computes the Wigner-3j symbol for given l1,l2,l3,m1,m2,m3. We
 * explicitly enforce the selection rules. */
inline double wigner3j(double l1, double l2, double l3, 
					 	double m1, double m2, double m3)
{
	// We enforce the selection rules.
	bool select(true);
	select = ( 
		   m1+m2+m3==0.0 
		&& std::floor(l1+l2+l3)==(l1+l2+l3)
		&& l3 >= std::fabs(l1-l2) 
		&& l3 <= l1+l2
		&& std::fabs(m1) <= l1
		&& std::fabs(m2) <= l2
		&& std::fabs(m3) <= l3
		);

	if (!select) return 0.0;

	// We compute the size of the resulting array.
	int size = (int)std::ceil(l2+l3-std::max(std::fabs(l2-l3),std::fabs(m1)))+1;

	// We prepare the output values.
	double l1min, l1max;
	double thrcof [size];
	int ierr;

	// External function call. 
	drc3jj_wrap(l2,l3,m2,m3,&l1min,&l1max,thrcof,size,&ierr);

	// We fetch and return the value with the proper l1 value.
	int index = (int)l1-l1min; 
	return thrcof[index];
}

/*! Computes the Clebsch-Gordan coefficient by relating it to the
 * Wigner 3j symbol. It sometimes eases the notation to use the 
 * Clebsch-Gordan coefficients directly. */
inline double clebschGordan(double l1, double l2, double l3,
							double m1, double m2, double m3)
{
	// We simply compute it via the 3j symbol. 
	return (pow(-1.0,l1-l2+m3)*sqrt(2.0*l3+1.0)*wigner3j(l1,l2,l3,m1,m2,-m3));
}

/*! Computes the Wigner-6j symbol for given, l1, l2, l3, l4, l5, l6.
 * We explicitly enforce the selection rules. */
inline double wigner6j(double l1, double l2, double l3, 
						double l4, double l5, double l6)
{
	// We enforce the selection rules. 
	bool select(true);

	// Triangle relations for the four tryads
	select = ( 
		   std::fabs(l1-l2) <= l3 && l3 <= l1+l2
		&& std::fabs(l1-l5) <= l6 && l6 <= l1+l5
		&& std::fabs(l4-l2) <= l6 && l6 <= l4+l2
		&& std::fabs(l4-l5) <= l3 && l3 <= l4+l5
		);

	// Sum rule of the tryads
	select = (
		   std::floor(l1+l2+l3)==(l1+l2+l3)
		&& std::floor(l1+l5+l6)==(l1+l5+l6)
		&& std::floor(l4+l2+l6)==(l4+l2+l6)
		&& std::floor(l4+l5+l3)==(l4+l5+l3)
		);

	if (!select) return 0.0;

	// We compute the size of the resulting array. 
	int size = (int)std::ceil(std::min(l2+l3,l5+l6)-std::max(std::fabs(l2-l3),std::fabs(l5-l6)))+1;

	// We prepare the output values
	double l1min, l1max;
	double sixcof [size];
	int ierr;

	// External function call
	drc6j_wrap(l2,l3,l4,l5,l6,&l1min,&l1max,sixcof,size,&ierr);

	// We fetch and return the coefficient with the proper l1 value. 
	int index = (int)l1-l1min;
	return sixcof[index];
}

template <typename T> double sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
}

#endif  // WIGNER_SYMBOLS_H