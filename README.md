Wigner Symbols
==============
[![Build status](https://travis-ci.org/valandil/wignerSymbols.svg?branch=master)](https://travis-ci.org/valandil/wignerSymbols)
[![Code coverage report](https://codecov.io/gh/valandil/wignerSymbols/branch/master/graph/badge.svg)](https://codecov.io/gh/valandil/wignerSymbols)
[![Join the chat at https://gitter.im/valandil/wignerSymbols](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/valandil/wignerSymbols?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/5354/valandil/wignerSymbols.svg)](https://zenodo.org/badge/latestdoi/5354/valandil/wignerSymbols)

A C++ ensemble of functions to compute the Wigner 3j- and 6j- symbols. It reimplements the algorithm designed
by Schulten and Gordon in C++, but also contains the original Fortran implementation. 
It can act either as a complete C++ replacement to the original Fortran implementation, 
or a C++ interface to them. See the API docs for details. 
It can either compute an array of Wigner 3j or 6j symbols, or a single
coefficient. It also computes the Clebsch-Gordan coefficients.

Note that there is a third party [Python wrapper](https://github.com/jeffzhen/wignerpy) 
available on GitHub.


## Compilation Instructions
This library uses CMake to help the build process. First, download the source code. 
It is recommended to create a separate directory for building, i.e.
```bash 
mkdir build/
cd build/
```
Then, run
```bash
cmake .. && make && sudo make install
```
By default, the library is installed to `/usr/lib/` and the include files are in `/usr/include/`.
To install to another directory, say `/usr/local/`, use the command-line argument
```bash
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr/local && make && sudo make install
```

## API documentation
We list the user-facing functions that compute the Wigner symbols. The functions are
behind the namespace `WignerSymbols`. Both the C++ implementation (my own) and the
original Fortran implementation (by Schulten and Gordon) are provided by this package.
In my own tests, I had found that the original Fortran implementation provided only 
`single` precision. This might be compiler-dependent, so you might have better luck. 
The Fortran version is also somewhat faster (10% faster, approximately). In any case, 
this program provides either a complete C++ replacement, or a C++ interface to the 
evaluation of Wigner symbols. 

### C++ implementation

  + `std::vector<double> wigner3j(double l2, double l3, double m1, double m2, double m3)`<br />
    Computes Wigner 3j symbols with all possible values of `l1`. Returns an `std::vector<double>` with the 
    coefficients sorted by increasing values of `l1`.
  + `double wigner3j(double l1, double l2, double l3, double m1, double m2, double m3)`<br />
    Computes a specific Wigner 3j symbol. 
  + `double clebschGordan(double l1, double l2, double l3, double m1, double m2, double m3)`<br />
    Computes a specific Clebch-Gordan coeffcient.
  + `std::vector<double> wigner6j(double l2, double l3, double l4, double l5, double l6)`<br />
    Computes Wigner 6j symbols with all possible values of `l1`. Returns an `std::vector<double>` with the 
    coefficients sorted by increasing values of `l1`.
  + `double wigner6j(double l1, double l2, double l3, double l4, double l5, double l6)`<br />
    Computes a specific Wigner 6j symbol.

### Fortran implementation

  + `std::vector<double> wigner3j_f(double l2, double l3, double m1, double m2, double m3)` <br />
    Computes Wigner 3j symbols with all possible values of `l1`. 
  + `double wigner3j_f(double l1, double l2, double l3, double m1, double m2, double m3)`<br />
    Computes a specific Wigner 3j symbol. 
  + `double clebschGordan_f(double l1, double l2, double l3, double m1, double m2, double m3)`<br />
    Computes a specific Clebch-Gordan coeffcient.
  + `std::vector<double> wigner6j_f(double l2, double l3, double l4, double l5, double l6)` <br />
    Computes Wigner 6j symbols with all possible values of `l1`. 
  + `double wigner6j_f(double l1, double l2, double l3, double l4, double l5, double l6)`<br />
    Computes a specific Wigner 6j symbol.

## Bibliography 
  + K. Schulten and R. G. Gordon, _Recursive evaluation of 3j and 6j coefficients_, Comput. Phys. Commun. **11**, 269–278 (1976). DOI: [10.1016/0010-4655(76)90058-8](https://dx.doi.org/10.1016/0010-4655(76)90058-8)
  + K. Schulten, _Exact recursive evaluation of 3j- and 6j-coefficients for quantum-mechanical coupling of angular momenta_, J. Math. Phys. **16**, 1961 (1975). DOI: [10.1063/1.522426](https://dx.doi.org/10.1063/1.522426).
  + J. Luscombe and M. Luban, _Simplified recursive algorithm for Wigner 3j and 6j symbols_, Phys. Rev. E **57**, 7274–7277 (1998). DOI: [10.1103/PhysRevE.57.7274](https://dx.doi.org/10.1103/PhysRevE.57.7274).

