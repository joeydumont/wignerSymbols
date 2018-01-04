/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the LGPL was not -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/lgpl.html.              -/
 ********************************************************/

#include "../include/wignerSymbols/wignerSymbols-fortran.h"

namespace WignerSymbols {


/*! Computes a string of Wigner-3j symbols for given l2, l3, m1, m2, m3. */
std::vector<double> wigner3j_f(double l2, double l3, double m1, double m2, double m3)
{
  // We prepare the size of the resulting array.
  int size = (int)std::ceil(l2+l3-std::max(std::fabs(l2-l3),std::fabs(m1)))+1;

  // We prepare the output values.
  double l1min, l1max;
  std::vector<double> thrcof(size);
  int ierr;

  // External function call.
  drc3jj_wrap(l2,l3,m2,m3,&l1min,&l1max,thrcof.data(),size,&ierr);

  return thrcof;
}

/*! Computes the Wigner-3j symbol for given l1,l2,l3,m1,m2,m3. We
 * explicitly enforce the selection rules. */
double wigner3j_f(double l1, double l2, double l3,
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
  std::vector<double> thrcof(size);
  int ierr;

  // External function call.
  drc3jj_wrap(l2,l3,m2,m3,&l1min,&l1max,thrcof.data(),size,&ierr);

  // We fetch and return the value with the proper l1 value.
  int index = (int)(l1-l1min);
  return thrcof[index];
}

double wigner6j_f(double l1, double l2, double l3,
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
  std::vector<double> sixcof(size);
  int ierr;

  // External function call
  drc6j_wrap(l2,l3,l4,l5,l6,&l1min,&l1max,sixcof.data(),size,&ierr);

  // We fetch and return the coefficient with the proper l1 value.
  int index = (int)(l1-l1min);
  return sixcof[index];
}
}
