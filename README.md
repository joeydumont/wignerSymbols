Wigner Symbols
==============

A C++ ensemble of functions to compute the Wigner 3j- and 6j- symbols. It implements the algorihtm designed
by Schulten and Gordon. It can either compute an array of Wigner 3j or 6j symbols, or a single
coefficient. It also computes the Clebsch-Gordan coefficients.

## API documentation
We list the user-facing functions that compute the Wigner symbols. The functions are
behind the namespace `WignerSymbols`.

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


## Bibliography 
  + K. Schulten and R. G. Gordon, _Recursive evaluation of 3j and 6j coefficients_, Comput. Phys. Commun. **11**, 269–278 (1976).
  + K. Schulten, _Exact recursive evaluation of 3j- and 6j-coefficients for quantum-mechanical coupling of angular momenta,_ J. Math. Phys. **16**, 1961 (1975).

