/*******************************************************-/
 * This source code is subject to the terms of the GNU  -/
 * Lesser Public License. If a copy of the LGPL was not -/
 * distributed with this file, you can obtain one at    -/
 * https://www.gnu.org/licenses/lgpl.html.              -/
 ********************************************************/

/** \file commonFunctions.h
 *
 *  \author Joey Dumont <joey.dumont@gmail.com>
 *
 *  \since 2014-07-10
 *
 *  \brief Defines some common functions to the C++ and Fortran code.
 *
 * We define some functions that will be of use in both the C++ and Fortran
 * parts of this library.
 *
 */

namespace WignerSymbols {

template <typename T>
double sgn(T val)
{
    int sgn = (T(0) < val) - (val < T(0));
    if (sgn == 0)
        return 1.0;
    else
        return (double)sgn;
}

} // namespace WignerSymbols
