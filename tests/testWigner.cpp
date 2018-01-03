/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the	LGPL was not-/
 * distributed with this file, you	can obtain one at   -/
 * https://www.gnu.org/licenses/lgpl.html				-/
 ********************************************************/ 

#include "wignerSymbols.h"

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <time.h>

#define EPS 1.0e-5

using namespace std;
using namespace arma;

// Orthogonality relations and sum rules for the 3j symbols [Brink & Satchler, pp. 136+139].
double firstOrthoRelation3j(double l1, double l2, double l3, double m3);
double sumOverM23j(double l1, double l2, double l3, double m3);

// Orthogonality relations for the 6j symbols [Brink & Satchler pp. 142-3].
double firstSumOverL3(double l1, double l2, double l6);
double seconSumOverL3(double l1, double l2, double l6);

// Timing computation of Wigner symbols. 
double timingWignerSymbols(double l2, double l3, double m1, double m2, double m3);
double timingWignerSymbolsF(double l2, double l3, double m1, double m2, double m3);

int main(int argc, char* argv[])
{
  // We evalaute ThreeJ((1,0),(2,0),(1,0))
  std::cout << WignerSymbols::wigner3j(1,2,1,0,0,0) << std::endl;
	// We time the computation for different sizes.
	int N = 500;
	double l2(500.0),m1(-10.0),m2(60.0),m3(-50.0);
	mat times(N,2);
	mat timesF(N,2);

	for (int i=0;i<N;i++)
	{
		int l3 = i+std::fabs(m3);
		times(i,0) = (int)std::ceil(l2+l3-std::max(std::fabs(l2-l3),std::fabs(m1)))+1;
		timesF(i,0) = times(i,0);
		times(i,1) = timingWignerSymbols(l2,l3,m1,m2,m3);
		timesF(i,1) = timingWignerSymbolsF(l2,l3,m1,m2,m3);
	}

	times.save("times.dat", raw_ascii);
	timesF.save("timesF.dat",raw_ascii);

	// Generate values for l1, l2 and derive the rest of the allowed values.
        if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <lMax>\n";
            return 1;
        }
	double lMax = atof(argv[1]);
	colvec lVec = linspace(0,lMax,lMax+1);

	// Run over l1
	for (int i=0;i<lMax+1;i++)
	{
		double l1 = lVec(i);

		// Run over l2
		for (int j=0;j<lMax+1;j++)
		{
			
			double l2 = lVec(j);
			double l3min = std::fabs(l1-l2);
			double l3max = l1+l2;

			// Run over l3
			for (int k=0;k<(l3max-l3min+1);k++)
			{
				double l3 = l3min+k;

				double test = firstSumOverL3(l1,l2,l3); 

				if (!test)
					cout << "The first sum over L3 test failed for \n\tl1 = " << l1 << "\n\tl2 = " << l2 << "\n\tl6 = " << l3 << endl;

				// Run over m3
				mat precisionWigner(2*l3+1,2);
				for (int l=0;l<(2*l3+1);l++)
				{
					double m3 = -l3+l;
					test = firstOrthoRelation3j(l1,l2,l3,m3);

					precisionWigner(l,0) = m3;
					precisionWigner(l,1) = test;

					if (!test)
						cout << "The first orthogonality test failed for \n\tl1 = " << l1 << "\n\tl2 = " << l2 << "\n\tl3 = " << l3 << "\n\tm3 = " << m3 << endl;

					if (l3!=0.0) test = sumOverM23j(l1,l2,l3,m3);

					if (!test)
						cout << "The sum over M2 test failed for \n\tl1 = " << l1 << "\n\tl2 = " << l2 << "\n\tl3 = " << l3 << "\n\tm3 = " << m3 << endl;
				}
				precisionWigner.save("precision.dat", raw_ascii);
			}

			double l6max = std::min(2.*l1,2.*l2);
			
			// Run over l6
			for (int k=0;k<l6max+1;k++)
			{
				bool test = seconSumOverL3(l1,l2,k);

				if (!test)
					cout << "The second sum over L3 test has failed for \n\tl1 = " << l1 << "\n\tl2 = " << l2 << "\n\tl6 = " << k << endl;
			}
		}
	}

		// We print some coefficients.
	std::vector<double> wigner6j = WignerSymbols::wigner6j(16,14,13,15,15);
	for (std::vector<double>::iterator it = wigner6j.begin(); it != wigner6j.end(); ++it)
	{
		std::cout << std::setprecision(10) << std::scientific <<  *it << std::endl;
	}
	std::cout << WignerSymbols::wigner6j(450,500,510,520,510,510) << std::endl;
	std::cout << WignerSymbols::wigner6j_f(450,500,510,520,510,510) << std::endl;

	return 0;
}

// Sum_{m1,m2} (2*l3+1)*Wigner3j(l1,l2,l3,m1,m2,m3)*3j(l1,l2,l3',m1,m2,m3') = delta(l3,l3')*delta(m3,m3')
double firstOrthoRelation3j(double l1, double l2, double l3, double m3)
{
	// We determine all possible values of m1 and m2.
	int m1Size = ((int)(2.0*l1+1.0));
	int m2Size = ((int)(2.0*l2+1.0)); 

	double sum = 0.0;
	for (int i=0;i<m1Size;i++)
	{
		for (int j=0;j<m2Size;j++)
		{
			sum += (2.0*l3+1.0)*WignerSymbols::wigner3j_f(l1,l2,l3,((double)(-l1+i)),((double)(-l2+j)),m3)
				*WignerSymbols::wigner3j_f(l1,l2,l3,((double)(-l1+i)),((double)(-l2+j)),m3);
		}
	}

	double diff = std::fabs(1.0-sum);

	return diff;
}

// sum_{m1,m2} m2*(2*l3+1)c_wigner3j(l1,l2,l3,m1,m2,-m3)^2 = m3*(l3*(l3+1)+l2*(l2+1)-l1*(l1+1))/(2*l3*(l3+1))
double sumOverM23j(double l1, double l2, double l3, double m3)
{
	double sum = 0.0;
	// We sum over all allowable values of m1+m2+m3=0.
	for (int i=0;i<(2*l1+1);i++)
	{
		double m1 = -l1+i;

		for (int j=0;j<(2*l2+1);j++)
		{
			double m2 = -l2+j;

			if (std::fabs(m1+m2-m3)<1.0e-6)
				sum += m2*(2.0*l3+1.0)*pow(WignerSymbols::wigner3j(l1,l2,l3,m1,m2,-m3),2.0);
		}
	}

	double value = m3*(l3*(l3+1.)+l2*(l2+1.)-l1*(l1+1.))/(2.*l3*(l3+1.));

	return std::fabs(sum-value);
}

// sum_{l3} = (-1)^(2*l3)*(2*l3+1)*Wigner6j(l1,l2,l3,l1,l2,l6) = 1
double firstSumOverL3(double l1, double l2, double l6)
{
	// We sum over the value of l3. 
	double l3min = std::fabs(l2-l1);
	double l3max = l1+l2;

	double sum = 0.0;
	for (int i=0;i<l3max-l3min+1;i++)
	{
		double l3 = l3min+i;
		sum += pow(-1.,2.0*l3)*(2.0*l3+1.0)*WignerSymbols::wigner6j(l1,l2,l3,l1,l2,l6);
	}

	double diff = std::fabs(1.0-sum);

	return diff;
}

// sum_{l3} = (-1)^(l1+l2+l3)*(2*l3+1)*Wigner6j(l1,l2,l3,l2,l1,l6) = delta(l6,0)*sqrt((2*l1+1)*(2*l2+1))
double seconSumOverL3(double l1, double l2, double l6)
{
		// We sum over the value of l3. 
	double l3min = std::fabs(l1-l2);
	double l3max = l1+l2;

	double sum = 0.0;
	for (int i=0;i<l3max-l3min+1;i++)
	{
		double l3 = l3min+i;
		sum += pow(-1.,l1+l2+l3)*(2.0*l3+1.0)*WignerSymbols::wigner6j(l1,l2,l3,l2,l1,l6);
	}

	double value = (l6==0.0 ? sqrt((2.*l1+1.)*(2.*l2+1.)) : 0.0) ;

	//double diff = std::fabs(value-sum);

	return std::fabs(value-sum);
	
}

double timingWignerSymbols(double l2, double l3, double m1, double m2, double m3)
{
	// We start the timer. 
	clock_t start = clock();

	// We call the subroutine.
	colvec test = WignerSymbols::wigner3j(l2,l3,m1,m2,m3);

	clock_t end = clock();

	return ((double)end-start)/CLOCKS_PER_SEC;
}

double timingWignerSymbolsF(double l2, double l3, double m1, double m2, double m3)
{
	// We start the timer. 
	clock_t start = clock();

	// We call the subroutine.
	WignerSymbols::wigner3j(l2+l3,l2,l3,m1,m2,m3);

	clock_t end = clock();

	return ((double)end-start)/CLOCKS_PER_SEC;
}

