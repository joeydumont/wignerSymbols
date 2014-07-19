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
  double test14 = WignerSymbols::wigner3j_f(325, 999, 1221, 280, 899, -1179);
  double test15 = WignerSymbols::wigner3j_f(999, 1221, 325,  899, -1179, 280);
  double test16 = WignerSymbols::wigner3j_f(1221, 325, 999, -1179, 280, 899);

  double test21 = WignerSymbols::wigner3j(693, 896, 1371, 513, -838, 325);
  double test22 = WignerSymbols::wigner3j(896, 1371, 693, -838, 325, 513);
  double test23 = WignerSymbols::wigner3j(1371, 693, 896, 325, 513, -838);
  double test24 = WignerSymbols::wigner3j_f(693, 896, 1371, 513, -838, 325);
  double test25 = WignerSymbols::wigner3j_f(896, 1371, 693, -838, 325, 513);
  double test26 = WignerSymbols::wigner3j_f(1371, 693, 896, 325, 513, -838);

  double test31 = WignerSymbols::wigner3j(772, 874, 1231, 442, 756, -1198);
  double test32 = WignerSymbols::wigner3j(874, 1231, 772, 756, -1198, 442);
  double test33 = WignerSymbols::wigner3j(1231, 772, 874, -1198, 442, 756);
  double test34 = WignerSymbols::wigner3j_f(772, 874, 1231, 442, 756, -1198);
  double test35 = WignerSymbols::wigner3j_f(874, 1231, 772, 756, -1198, 442);
  double test36 = WignerSymbols::wigner3j_f(1231, 772, 874, -1198, 442, 756);

  std::cout << test11 << std::endl;
  std::cout << WignerSymbols::wigner3j(693, 896, 1371, 513, -838, 325)  << std::endl;
  std::cout << WignerSymbols::wigner3j(772, 874, 1231, 442, 756, -1198) << std::endl;
  std::cout << WignerSymbols::wigner3j(842, 996, 1714, 376, 979, -1355) << std::endl;
  std::cout << WignerSymbols::wigner3j(101, 987, 1084, -80, -935, 1015) << std::endl;
  std::cout << WignerSymbols::wigner3j( 51, 1003, 978, -32,  993, -961) << std::endl;
  std::cout << WignerSymbols::wigner3j(217, 1008, 1107, -76, -987, 1063)<< std::endl;
  std::cout << WignerSymbols::wigner3j(408, 894, 954, -114, -840, 954)  << std::endl;
  std::cout << WignerSymbols::wigner3j(808, 980, 1734, 341, 954, -1295) << std::endl;
  std::cout << WignerSymbols::wigner3j(827, 1020, 1570, 590, 902, -1492)<< std::endl;

  std::cout << "C++ impl.: " << test11 << ", " << test12 << ", " << test13 << std::endl;
  std::cout << "FOR impl.: " << test14 << ", " << test15 << ", " << test16 << std::endl;
  std::cout << "C++ impl.: " << test21 << ", " << test22 << ", " << test23 << std::endl;
  std::cout << "FOR impl.: " << test24 << ", " << test25 << ", " << test26 << std::endl;
  std::cout << "C++ impl.: " << test31 << ", " << test32 << ", " << test33 << std::endl;
  std::cout << "FOR impl.: " << test34 << ", " << test35 << ", " << test36 << std::endl;

  std::vector<double> test4 = WignerSymbols::wigner3j(992, 1243, 196, -901, 705);
  double test41 = WignerSymbols::wigner3j(529, 992, 1243, 196, -901, 705);
  double test42 = WignerSymbols::wigner3j_f(529, 992, 1243, 196, -901, 705);
  std::cout << "C++ impl.: " << test41 << std::endl;
  std::cout << "FOR impl.: " << test42 << std::endl;
  std::cout << std::endl;

  return 0;
}
 
