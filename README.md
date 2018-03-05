# FPML
Fourth order Parallelizable Modified Laguerre method. This is a Fortran library collection of solvers and their dependencies for the roots of an univariate polynomial. 

## Authors
- [Thomas R. Cameron](https://thomasrcameron.com),
Davidson College, NC

## Instalation Instructions

## Tests
### Initial Estimate
Below is a graph of the initial estimates and exact roots of the polynomial <img src="https://latex.codecogs.com/svg.latex?\Large&space;p(z)=1+3000z+3000000z^{2}+1000000000z^{3}+z^{10}" title="\Large p(z)=1+3000z+3000000z^{2}+1000000000z^{3}+z^{10}" />.
![alt text](tests/figures/module.png?raw=true)

As can be seen from this graph, the initial estimates are incredibly close to the exact roots of the polynomial. While the outcome is not always this favorable, this example highlights the 
### Random Polynomials
Random complex polynomials whose coefficients are uniformly distributed over the interval [-1,1] are used to compare FPML against Polzeros and the singleshift version of AMVW. The plot below includes the elapsed time measured in seconds and the accuracy which is measured as the forward error. Iterations are run for polynomials of degree 100 to degree 6400, doubling the degree on each step. For each iteration there are 25 tests performed, the average time elapsed and the average of the maximum forward error over all these tests is recorded.
![alt text](tests/figures/rand_poly.png?raw=true)
### Roots of Unity
The polynomial <img src="https://latex.codecogs.com/svg.latex?\Large&space;z^{n}-1" title="\Large z^{n}-1" /> has as roots the n roots of unity: <img src="https://latex.codecogs.com/svg.latex?\Large&space;\cos(\frac{2\pi j}{n})" title="\Large \cos(\frac{2\pi j}{n})" /> for <img src="https://latex.codecogs.com/svg.latex?\Large&space;j=1,\ldots,n" title="\Large j=1,\ldots,n" />. These polynomials are used to compare FPML against Polzeros and the singleshift version of AMVW. The plot below includes the elapsed time measured in seconds and the accuracy which is measured as the absolute difference between the computed and exact roots. Iterations are run for polynomials of degree 100 to degree 6400, doubling the degree on each step. For each iteration there are 25 tests performed, the average time elapsed over the number of tests and the average error over the degree is recorded. 
![alt text](tests/figures/unity.png?raw=true)


## Related Articles
This software is based on the following articles:

1. O. Aberth, Iteration methods for finding all zeros of a polynomial simultaneously, Math. Comp. 27 (1973), no. 122, 339-344.