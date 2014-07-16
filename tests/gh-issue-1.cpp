/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the LGPL was not -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/lgpl.html               -/
 ********************************************************/

/*! \file gh-issue-1.cpp
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *  \since 2014-07-11
 *  \brief Tests Issue 1 of WignerSymbols (see valandil.github.com/wignerSymbols.
 *  \copyright LGPL
 * This file tests the bug described in Issue #1 of WignerSymbols (GitHub).
 * It seems that for some values of l1, the algorithm returns a set 
 * of zeros even though it should return proper coefficients. 
 * Permuting the l's lead to the proper output. 
 */

#include <wignerSymbols.h>

int main (int argc, char* argv[])
{
  double test11 = WignerSymbols::wigner3j(325, 999, 1221, 280, 899, -1179);
  double test12 = WignerSymbols::wigner3j(999, 1221, 325,  899, -1179, 280);
  double test13 = WignerSymbols::wigner3j(1221, 325, 999, -1179, 280, 899);

  double test21 = WignerSymbols::wigner3j(693, 896, 1371, 513, -838, 325);
  double test22 = WignerSymbols::wigner3j(896, 1371, 693, -838, 325, 513);
  double test23 = WignerSymbols::wigner3j(1371, 693, 896, 325, 513, -838);

  double test31 = WignerSymbols::wigner3j(772, 874, 1231, 442, 756, -1198);
  double test32 = WignerSymbols::wigner3j(874, 1231, 772, 756, -1198, 442);
  double test33 = WignerSymbols::wigner3j(1231, 772, 874, -1198, 442, 756);

//[ 772 874 1231 442 756 -1198] 0.00171146251356 0.0 0.0 0.00171146251356
//[ 842 996 1714 376 979 -1355] -9.12730017838e-14 -0.0 0.0 -9.12730017838e-14
//[ 101 987 1084 -80 -935 1015] -0.00274486614091 -0.0 0.0 -0.00274486614091
//[ 51 1003 978 -32 993 -961] -0.00561491486489 -0.0 0.0 -0.00561491486489
//[ 217 1008 1107 -76 -987 1063] -0.00129163704248 -0.0 0.0 -0.00129163704248
//[ 408 894 954 -114 -840 954] 0.000141186572123 0.0 0.0 0.000141186572123
//[ 808 980 1734 341 954 -1295] -1.13581014918e-35 -0.0 0.0 -1.13581014918e-35
//[ 827 1020 1570 590 902 -1492] -0.00091099434249 -0.0 0.0 -0.00091099434249

  std::cout << "C++ impl.: " << test11 << ", " << test12 << ", " << test13 << std::endl;
  std::cout << "C++ impl.: " << test21 << ", " << test22 << ", " << test23 << std::endl;
  std::cout << "C++ impl.: " << test31 << ", " << test32 << ", " << test33 << std::endl;

  std::cout << WignerSymbols::maxThree<double>(2.0, 1.0, 4.0) << std::endl;

  return 0;
}
 
