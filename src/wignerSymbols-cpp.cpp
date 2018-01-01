/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the LGPL was not -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/lgpl.html.              -/
 ********************************************************/

#include "../include/wignerSymbols/commonFunctions.h"
#include "../include/wignerSymbols/wignerSymbols-cpp.h"

namespace WignerSymbols {
std::vector<double> wigner3j(double l2, double l3,
			     double m1, double m2, double m3)
{
	// We compute the numeric limits of double precision.
	double huge = sqrt(std::numeric_limits<double>::max()/20.0);
	double srhuge = sqrt(huge);
	double tiny = std::numeric_limits<double>::min();
	double srtiny = sqrt(tiny);
	double eps = std::numeric_limits<double>::epsilon();

	// We enforce the selection rules.
	bool select(true);
	select = (
		   std::fabs(m1+m2+m3)<eps
		&& std::fabs(m2) <= l2+eps
		&& std::fabs(m3) <= l3+eps
		);

	if (!select) return std::vector<double>(1,0.0);

	// We compute the limits of l1.
	double l1min = std::max(std::fabs(l2-l3),std::fabs(m1));
	double l1max = l2+l3;

	// We compute the size of the resulting array.
	int size = (int)std::floor(l1max-l1min+1.0+eps);
	std::vector<double> thrcof(size,0.0);

	// If l1min=l1max, we have an analytical formula.
	if (size==1)
	{
		thrcof[0] = pow(-1.0,std::floor(std::fabs(l2+m2-l3+m3)))/sqrt(l1min+l2+l3+1.0);
	}

	// Another special case where the recursion relation fails.
	else
	{
		// We start with an arbitrary value.
		thrcof[0] = srtiny;

		// From now on, we check the variation of |alpha(l1)|.
		double alphaNew, l1(l1min);
		if (l1min==0.0)
			alphaNew = -(m3-m2+2.0*wigner3j_auxB(l1,l2,l3,m1,m2,m3))/wigner3j_auxA(1.0,l2,l3,m1,m2,m3);
		else
			alphaNew = -wigner3j_auxB(l1min,l2,l3,m1,m2,m3)
							/(l1min*wigner3j_auxA(l1min+1.0,l2,l3,m1,m2,m3));

		// We compute the two-term recursion.
		thrcof[1] = alphaNew*thrcof[0];

		// If size > 2, we start the forward recursion.
		if (size>2)
		{
			// We start with an arbitrary value.
			thrcof[0] = srtiny;

			// From now on, we check the variation of |alpha(l1)|.
			double alphaOld, alphaNew, beta, l1(l1min);
			if (l1min==0.0)
				alphaNew = -(m3-m2+2.0*wigner3j_auxB(l1,l2,l3,m1,m2,m3))/wigner3j_auxA(1.0,l2,l3,m1,m2,m3);
			else
				alphaNew = -wigner3j_auxB(l1min,l2,l3,m1,m2,m3)
							/(l1min*wigner3j_auxA(l1min+1.0,l2,l3,m1,m2,m3));

			// We compute the two-term recursion.
			thrcof[1] = alphaNew*thrcof[0];

			// We compute the rest of the recursion.
			int i = 1;
			bool alphaVar = false;
			do
			{
				// Bookkeeping:
				i++;					// Next term in recursion
				alphaOld = alphaNew;	// Monitoring of |alpha(l1)|.
				l1 += 1.0;				// l1 = l1+1

				// New coefficients in recursion.
				alphaNew = -wigner3j_auxB(l1,l2,l3,m1,m2,m3)
							/(l1*wigner3j_auxA(l1+1.0,l2,l3,m1,m2,m3));

				beta = -(l1+1.0)*wigner3j_auxA(l1,l2,l3,m1,m2,m3)
						/(l1*wigner3j_auxA(l1+1.0,l2,l3,m1,m2,m3));

				// Application of the recursion.
				thrcof[i] = alphaNew*thrcof[i-1]+beta*thrcof[i-2];

				// We check if we are overflowing.
				if (std::fabs(thrcof[i])>srhuge)
				{
					std::cout << "We renormalized the forward recursion." << std::endl;
					for (std::vector<double>::iterator it = thrcof.begin(); it != thrcof.begin()+i; ++it)
					{
						//if (std::fabs(*it) < srtiny) *it = 0;
						//else
						*it /= srhuge;
					}
				}

				// This piece of code checks whether we have reached
				// the classical region. If we have, the second if
				// sets alphaVar to true and we break this loop at the
				// next iteration because we need thrcof(l1mid+1) to
				// compute the scalar lambda.
				if (alphaVar) break;

				if (std::fabs(alphaNew)-std::fabs(alphaOld)>0.0)
					alphaVar=true;

			}	while (i<(size-1));	// Loop stops when we have computed all values.

			// If this is the case, we have stumbled upon a classical region.
			// We start the backwards recursion.
			if (i!=size-1)
			{
				// We keep the two terms around l1mid to compute the factor later.
				double l1midm1(thrcof[i-2]),l1mid(thrcof[i-1]),l1midp1(thrcof[i]);

				// We compute the backward recursion by providing an arbitrary
				// startint value.
				thrcof[size-1] = srtiny;

				// We compute the two-term recursion.
				l1 = l1max;
				alphaNew = -wigner3j_auxB(l1,l2,l3,m1,m2,m3)
								/((l1+1.0)*wigner3j_auxA(l1,l2,l3,m1,m2,m3));
				thrcof[size-2] = alphaNew*thrcof[size-1];

				// We compute the rest of the backward recursion.
				int j = size-2;
				do
				{
					// Bookkeeping
					j--;			// Previous term in recursion.
					l1 -= 1.0;		// l1 = l1-1

					// New coefficients in recursion.
					alphaNew = -wigner3j_auxB(l1,l2,l3,m1,m2,m3)
								/((l1+1.0)*wigner3j_auxA(l1,l2,l3,m1,m2,m3));
					beta = -l1*wigner3j_auxA(l1+1.0,l2,l3,m1,m2,m3)
									/((l1+1.0)*wigner3j_auxA(l1,l2,l3,m1,m2,m3));

					// Application of the recursion.
					thrcof[j] = alphaNew*thrcof[j+1]+beta*thrcof[j+2];

					// We check if we are overflowing.
					if (std::fabs(thrcof[j]>srhuge))
					{
						std::cout << "We renormalized the backward recursion." << std::endl;
						for (std::vector<double>::iterator it = thrcof.begin()+j; it != thrcof.end(); ++it)
						{
							//if (std::fabs(*it) < srtiny) *it = 0;
							//else
							*it /= srhuge;
						}
					}

				} while (j>(i-2)); // Loop stops when we are at l1=l1mid-1.

				// We now compute the scaling factor for the forward recursion.
				double lambda = (l1midp1*thrcof[j+2]+l1mid*thrcof[j+1]+l1midm1*thrcof[j])
										/(l1midp1*l1midp1+l1mid*l1mid+l1midm1*l1midm1);

				// We scale the forward recursion.
				for (std::vector<double>::iterator it = thrcof.begin(); it != thrcof.begin()+j; ++it)
				{
					*it *= lambda;
				}
			}
		}
	}

	// We compute the overall factor.
	double sum = 0.0;
	for (int k=0;k<size;k++)
	{
		sum += (2.0*(l1min+k)+1.0)*thrcof[k]*thrcof[k];
	}
	//std::cout << sum << std::endl;

	//std::cout << "(-1)^(l2-l3-m1): " << pow(-1.0,l2-l3-m1) << " sgn:" << sgn(thrcof[size-1]) << std::endl;
	double c1 = pow(-1.0,l2-l3-m1)*sgn(thrcof[size-1]);
	//std::cout << "c1: " << c1 << std::endl;
	for (std::vector<double>::iterator it = thrcof.begin(); it != thrcof.end(); ++it)
	{
		//std::cout << *it << ", " << c1 << ", ";
		*it *= c1/sqrt(sum);
		//std::cout << *it << std::endl;
	}
	return thrcof;
}

double wigner3j(double l1, double l2, double l3,
					double m1, double m2, double m3)
{
	// We enforce the selection rules.
	bool select(true);
	select = (
		   std::fabs(m1+m2+m3)<1.0e-10
		&& std::floor(l1+l2+l3)==(l1+l2+l3)
		&& l3 >= std::fabs(l1-l2)
		&& l3 <= l1+l2
		&& std::fabs(m1) <= l1
		&& std::fabs(m2) <= l2
		&& std::fabs(m3) <= l3
		);

	if (!select) return 0.0;

	// We compute l1min and the position of the array we will want.
	double l1min = std::max(std::fabs(l2-l3),std::fabs(m1));

	// We fetch the proper value in the array.
	int index = (int)(l1-l1min);

	return wigner3j(l2,l3,m1,m2,m3)[index];
}

std::vector<double> wigner6j(double l2, double l3,
					double l4, double l5, double l6)
{
	// We compute the numeric limits of double precision.
	double huge = std::numeric_limits<double>::max();
	double srhuge = sqrt(huge);
	double tiny = std::numeric_limits<double>::min();
	double srtiny = sqrt(tiny);
	double eps = std::numeric_limits<double>::epsilon();

	// We enforce the selection rules.
	bool select(true);

	// Triangle relations for the four tryads
	select = (
		std::fabs(l4-l2) <= l6 && l6 <= l4+l2
		&& std::fabs(l4-l5) <= l3 && l3 <= l4+l5
		);

	// Sum rule of the tryads
	select = (
		std::floor(l4+l2+l6)==(l4+l2+l6)
		&& std::floor(l4+l5+l3)==(l4+l5+l3)
		);

	if (!select) return std::vector<double>(1,0.0);

	// We compute the limits of l1.
	double l1min = std::max(std::fabs(l2-l3),std::fabs(l5-l6));
	double l1max = std::min(l2+l3,l5+l6);

	// We compute the size of the resulting array.
	unsigned int size = (int)std::floor(l1max-l1min+1.0+eps);
	std::vector<double> sixcof(size,0.0);

	// If l1min=l1max, we have an analytical formula.
	if (size==1)
	{
		sixcof[0] = 1.0/sqrt((l1min+l1min+1.0)*(l4+l4+1.0));
		sixcof[0] *= ((int)std::floor(l2+l3+l5+l6+eps) & 1 ? -1.0 : 1.0);
	}

	// Otherwise, we start the forward recursion.
	else
	{
		// We start with an arbitrary value.
		sixcof[0] = srtiny;

		// From now on, we check the variation of |alpha(l1)|.
		double alphaNew, l1(l1min);

		if (l1min==0)
			alphaNew = -(l2*(l2+1.0)+l3*(l3+1.0)+l5*(l5+1.0)+l6*(l6+1.0)-2.0*l4*(l4+1.0))/wigner6j_auxA(1.0,l2,l3,l4,l5,l6);

		else
			alphaNew = -wigner6j_auxB(l1,l2,l3,l4,l5,l6)
							/(l1min*wigner6j_auxA(l1+1.0,l2,l3,l4,l5,l6));

		// We compute the two-term recursion.
		sixcof[1] = alphaNew*sixcof[0];

		if (size>2)
		{
			// We start with an arbitrary value.
			sixcof[0] = srtiny;

			// From now on, we check the variation of |alpha(l1)|.
			double alphaOld, alphaNew, beta, l1(l1min);
			if (l1min==0)
				alphaNew = -(l2*(l2+1.0)+l3*(l3+1.0)+l5*(l5+1.0)+l6*(l6+1.0)-2.0*l4*(l4+1.0))/wigner6j_auxA(1.0,l2,l3,l4,l5,l6);

			else
				alphaNew = -wigner6j_auxB(l1,l2,l3,l4,l5,l6)
							/(l1min*wigner6j_auxA(l1+1.0,l2,l3,l4,l5,l6));

			// We compute the two-term recursion.
			sixcof[1] = alphaNew*sixcof[0];

			// We compute the rest of the recursion.
			unsigned int i = 1;
			bool alphaVar = false;
			do
			{
				// Bookkeeping:
				i++;					// Next term in recursion
				alphaOld = alphaNew;	// Monitoring of |alpha(l1)|.
				l1 += 1.0;				// l1 = l1+1

				// New coefficients in recursion.
				alphaNew = -wigner6j_auxB(l1,l2,l3,l4,l5,l6)
							/(l1*wigner6j_auxA(l1+1.0,l2,l3,l4,l5,l6));

				beta = -(l1+1.0)*wigner6j_auxA(l1,l2,l3,l4,l5,l6)
						/(l1*wigner6j_auxA(l1+1.0,l2,l3,l4,l5,l6));

				// Application of the recursion.
				sixcof[i] = alphaNew*sixcof[i-1]+beta*sixcof[i-2];

				// We check if we are overflowing.
				if (std::fabs(sixcof[i]>srhuge))
				{
					std::cout << "We renormalized the forward recursion." << std::endl;
					for (std::vector<double>::iterator it = sixcof.begin(); it != sixcof.begin()+i; ++it)
					{
						*it /= srhuge;
					}
				}

				// This piece of code checks whether we have reached
				// the classical region. If we have, the second if
				// sets alphaVar to true and we break this loop at the
				// next iteration because we need sixcof(l1mid+1) to
				// compute the scalar.
				if (alphaVar) break;

				if (std::fabs(alphaNew)-std::fabs(alphaOld)>0.0)
					alphaVar=true;

			}	while (i<(size-1));	// Loop stops when we have computed all values.

			// If this is the case, we have stumbled upon a classical region.
			// We start the backwards recursion.
			if (i!=size-1)
			{
				// We keep the two terms around l1mid to compute the factor later.
				double l1midm1(sixcof[i-2]),l1mid(sixcof[i-1]),l1midp1(sixcof[i]);

				// We compute the backward recursion by providing an arbitrary
				// startint value.
				sixcof[size-1] = srtiny;

				// We compute the two-term recursion.
				l1 = l1max;
				alphaNew = -wigner6j_auxB(l1,l2,l3,l4,l5,l6)
								/((l1+1.0)*wigner6j_auxA(l1,l2,l3,l4,l5,l6));
				sixcof[size-2] = alphaNew*sixcof[size-1];

				// We compute the rest of the backward recursion.
				unsigned int j = size-2;
				do
				{
					// Bookkeeping
					j--;			// Previous term in recursion.
					l1 -= 1.0;		// l1 = l1-1

					// New coefficients in recursion.
					alphaNew = -wigner6j_auxB(l1,l2,l3,l4,l5,l6)
								/((l1+1.0)*wigner6j_auxA(l1,l2,l3,l4,l5,l6));
					beta = -l1*wigner6j_auxA(l1+1.0,l2,l3,l4,l5,l6)
								/((l1+1.0)*wigner6j_auxA(l1,l2,l3,l4,l5,l6));

					// Application of the recursion.
					sixcof[j] = alphaNew*sixcof[j+1]+beta*sixcof[j+2];

					// We check if we are overflowing.
					if (std::fabs(sixcof[j]>srhuge))
					{
						std::cout << "We renormalized the backward recursion." << std::endl;
						for (std::vector<double>::iterator it = sixcof.begin()+j; it != sixcof.end(); ++it)
						{
							*it /= srhuge;
						}
					}

				} while (j>(i-2)); // Loop stops when we are at l1=l1mid-1.

				// We now compute the scaling factor for the forward recursion.
				double lambda = (l1midp1*sixcof[j+2]+l1mid*sixcof[j+1]+l1midm1*sixcof[j])
									/(l1midp1*l1midp1+l1mid*l1mid+l1midm1*l1midm1);

				// We scale the forward recursion.
				for (std::vector<double>::iterator it = sixcof.begin(); it != sixcof.begin()+j; ++it)
				{
					*it *= lambda;
				}
			}
		}
	}

	// We compute the overall factor.
	double sum = 0.0;
	for (unsigned int k=0;k<size;k++)
	{
		sum += (2.0*(l1min+k)+1.0)*(2.0*l4+1.0)*sixcof[k]*sixcof[k];
	}
	double c1 = pow(-1.0,std::floor(l2+l3+l5+l6+eps))*sgn(sixcof[size-1])/sqrt(sum);

	for (std::vector<double>::iterator it = sixcof.begin(); it != sixcof.end(); ++it)
	{
		*it *= c1;
	}
	return sixcof;
}

double wigner6j(double l1, double l2, double l3,
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

	// We compute l1min and the position of the array we will want.
	double l1min = std::max(std::fabs(l2-l3),std::fabs(l5-l6));
	int index = (int)(l1-l1min);

	return wigner6j(l2,l3,l4,l5,l6)[index];
}

double wigner3j_auxA(double l1, double l2, double l3,
                                               double m1, double /*m2*/, double /*m3*/)
{
	double T1 = l1*l1-pow(l2-l3,2.0);
	double T2 = pow(l2+l3+1.0,2.0)-l1*l1;
	double T3 = l1*l1-m1*m1;

	return sqrt(T1*T2*T3);
}

double wigner3j_auxB(double l1, double l2, double l3,
						double m1, double m2, double m3)
{
	double T1 = -(2.0*l1+1.0);
	double T2 = l2*(l2+1.0)*m1;
	double T3 = l3*(l3+1.0)*m1;
	double T4 = l1*(l1+1.0)*(m3-m2);

	return T1*(T2-T3-T4);
}

double wigner6j_auxA(double l1, double l2, double l3,
						double /*l4*/, double l5, double l6)
{
	double T1 = l1*l1-pow(l2-l3,2.0);
	double T2 = pow(l2+l3+1.0,2.0)-l1*l1;
	double T3 = l1*l1-pow(l5-l6,2.0);
	double T4 = pow(l5+l6+1.0,2.0)-l1*l1;

	return sqrt(T1*T2*T3*T4);
}

double wigner6j_auxB(double l1, double l2, double l3,
						double l4, double l5, double l6)
{
	double T0 = 2.0*l1+1.0;

	double T1 = l1*(l1+1.0);
	double T2 = -l1*(l1+1.0)+l2*(l2+1.0)+l3*(l3+1.0);

	double T3 = l5*(l5+1.0);
	double T4 = l1*(l1+1.0)+l2*(l2+1.0)-l3*(l3+1.0);

	double T5 = l6*(l6+1.0);
	double T6 = l1*(l1+1.0)-l2*(l2+1.0)+l3*(l3+1.0);

	double T7 = 2.0*l1*(l1+1.0)*l4*(l4+1.0);

	return (T0*(T1*T2+T3*T4+T5*T6-T7));
}
}
