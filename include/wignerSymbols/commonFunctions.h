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
    return (T(0) < val) - (val < T(0));
}

template <typename T>
int maxThree(T a, T b, T c)
{
    int maxIndex;
    T max;
    if (a > b)
    {
        max = a;
        maxIndex = 1;
    }
    else 
    {
        max = b;
        maxIndex = 2;
    }

    if (max < c) 
    {
        max = c;
        maxIndex = 3;
    }

    return maxIndex;
}
}
