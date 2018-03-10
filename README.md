# FPML
Fourth order Parallelizable Modified Laguerre method. This is a Fortran library that includes a module of subroutines for computing the roots of an univariate polynomial. This module is name fpml and located in the src directory. 

## Authors
* [Thomas R. Cameron](https://thomasrcameron.com)
    * Mathematics and Computer Science Department, Davidson College
    * [Email: thcameron@davidson.edu](mailto:thcameron@davidson.edu)

## Installation Instructions
The main driver can be installed by running the install\_driver.sh script. When running the fpml\_driver, the program will expect the name of a file located in data\_files. This file is expected to store the coefficients of a polynomial as follows: The first row contains the degree and every subsequent line contains the coefficients from constant to leading (see poly1.dat for an example). If no data file is given, then the program will run on a random polynomial. The results including each computed root its backward error and condition number is recorded in the file results.dat.

All tests can be installed by running the install\_tests.sh script inside of the tests folder. If relevant, each test program has default values for the startDegree, endDegree, and maxit. If desired, these defaults can be overridden by providing these values when executing the program (see the run\_tests.sh script for an example on how to run each test executable).

After running a test the results will be in a file (named after the test) located in tests/data\_files. A figure displaying the data can be generated by compiling the associated TeX file (also named after the test). By running the script compile\_tex.sh all figures are generated and moved into the tests/figures directory.

We note that all installer scripts use *gfortran* to compile fortran files and *pdflatex* to compile tex files. When compiling fortran files the -O2 GNU optimization flag is used, and when compiling the pdflatex files the -shell-escape option is used for converting pdf files to png files, which requires *imagemagick* to be installed. Finally, all file paths used are relative to the directory structure of the GitHub directory FPML.
## Tests
The following numerical experiments are provided to motivate our choices in several design aspects of our algorithm, including our use of the monotone chain algorithm and our modification of Laguerre's method, and in addition to compare our algorithm with Polzeros [5] and AMVW [4]. All tests are run on an Intel Core i5 CPU running at 2.7 GHz with 16GB of memory.
### Initial Estimate
Below is a graph of the initial estimates and exact roots of the polynomial <img src="https://latex.codecogs.com/svg.latex?\Large&space;p(z)=1+3000z+3000000z^{2}+1000000000z^{3}+z^{10}" title="\Large p(z)=1+3000z+3000000z^{2}+1000000000z^{3}+z^{10}" />.

![alt text](tests/figures/init_est_acc.png?raw=true)

As can be seen from this graph, the initial estimates are incredibly close to the exact roots of the polynomial. While the outcome is not always this favorable, this example highlights the power of the method we use for computing initial estimates, originally due to Bini [5]. In addition, the plot below compares the divide and conquer method used by Bini and the monotone chain algorithm we employ, originally proposed in [3], for computing the upper envelope of the convex hull needed to compute the initial estimates.

![alt text](tests/figures/init_est_time.png?raw=true)

As expected, the monotone chain algorithm is significantly less expensive than the divide and conquer method for large degree polynomials. 
### Modifications
Random complex polynomials whose coefficients are uniformly distributed over the interval [-1,1] are used to compare several modifications of Laguerre's method and Newton's method. The modification of Newton's method was proposed in [1] and is referenced as Aberth's method. The modification we propose is referenced as concurrent Laguerre (Con. Lag.), and the modification used in [7,8] is referenced as sequential Laguerre (Seq. Lag.). The plot below includes the elapsed time measured in seconds and the accuracy which is measured as the maximum relative forward error. Iterations are run for polynomials of degree 80 to degree 10240, doubling the degree on each step. For each iteration there are 25 tests performed, the average time elapsed and the average accuracy over all these tests is recorded.

![alt text](tests/figures/methods.png?raw=true)

As can be seen from the graph above, all three modifications are comparable until the last two degrees tested. In these cases Seq. Lag. is failing to meet the stopping criterion for at least one root approximation, so the maximum backward error and condition number are much higher than the rest. In fact, we've observed that Seq. Lag. requires an increasing number of iterations as the degree of the polynomial grows. This, combined with Seq. Lag. not being parallelizable, is our motivation for using the Con. Lag. modification.
### Convergence
Three special polynomials are used to test the theoretical convergence properties of FPML. The first polynomial is <img src="https://latex.codecogs.com/svg.latex?\Large&space;z^{5}-1" title="\Large z^{5}-1" />, the second is the degree 10 Chebyshev polynomial, and the third polynomial is <img src="https://latex.codecogs.com/svg.latex?\Large&space;z^{20}+z^{19}+\cdots+z+1">. The error is measured as the maximum relative forward error. For each polynomial, the error after each iteration is recorded in the table below. The column Error-1 corresponds to the error in the roots approximations for the first polynomial, Error-2 for the second polynomial, and Error-3 for the third polynomial.

![alt text](tests/figures/conv.png?raw=true)

Note that evidence of fourth order convergence can be seen in each column. This is somewhat of a surprising feature of our method. Often, higher order methods do not display their convergence rate in practice. For instance, the method in [6] has fourth order convergence, but the radius of convergence is so small that by the time the estimates begin to converge the fourth order property will not be noticed in double precision floating arithmetic.
### Random Polynomials
Random complex polynomials whose coefficients are uniformly distributed over the interval [-1,1] are used to compare FPML against Polzeros and the singleshift version of AMVW. The plot below includes the elapsed time measured in seconds and the accuracy which is measured as the maximum relative forward error. Iterations are run for polynomials of degree 80 to degree 10240, doubling the degree on each step. For each iteration there are 25 tests performed, the average time elapsed and the average accuracy over all these tests is recorded.

![alt text](tests/figures/rand_poly.png?raw=true)

It is clear that FPML and Polzeros are very close with respect to time, with FPML having the slight advantage and AMVW being about 2 times slower on average. With respect to accuracy, FPML and AMVW are comparable, with AMVW have the slight advantage and Polzeros having 2 times larger accuracy on average.
### Roots of Unity
The polynomial <img src="https://latex.codecogs.com/svg.latex?\Large&space;z^{n}-1" title="\Large z^{n}-1" /> has as roots the n roots of unity: <img src="https://latex.codecogs.com/svg.latex?\Large&space;\cos(\frac{2\pi}{n}j)+i\sin(\frac{2\pi}{n}j)" title="\Large \cos(\frac{2\pi}{n}j)+i\sin(\frac{2\pi}{n}j)" /> for <img src="https://latex.codecogs.com/svg.latex?\Large&space;j=1,\ldots,n" title="\Large j=1,\ldots,n" />. These polynomials are used to compare FPML against Polzeros and the singleshift version of AMVW. The plot below includes the elapsed time measured in seconds and the accuracy which is measured as the average absolute difference between the computed and exact roots. Iterations are run for polynomials of degree 80 to degree 10240, doubling the degree on each step. For each iteration there are 25 tests performed, the average time elapsed over the number of tests and the accuracy is recorded. 

![alt text](tests/figures/unity.png?raw=true)

It is clear from the figure above that FPML and Polzeros are very close with respect to both time and accuracy. Note that the difference between FPML and Polzeros with respect to accuracy is within double precision unit roundoff and can therefore be attributed as noise. In these tests, AMVW is both slower and less accurate.
### Special Polynomials
Below is a table of the special polynomials used for providing additional comparisons between FPML, Polzeros, and the singleshift version of AMVW.

![alt text](tests/figures/spec_poly_list.png?raw=true)

For each special polynomial, the roots are computed using FPML, Polzeros, and AMVW. The maximum relative forward error is recorded in the table below. 

![alt text](tests/figures/spec_poly_results.png?raw=true)

It is clear that FPML and Polzeros have comparable accuracy when solving for the roots of each of the above special polynomials, often FPML is better by an order of magnitude. In addition, AMVW is performing much worse than both FPML and Polzeros in many of these tests.
## References

1. O. Aberth, *Iteration methods for finding all zeros of a polynomial simultaneously*, Math. Comp. 27 (1973), no. 122, 339-344.
2. F. S. Acton, *Numerical methods that work*, Harper and Row, New York, New York, 1970.
3. A. M. Andrew, *Another efficient algorithm for convex hulls in two dimensions*, Info. Proc. Letters 9 (1979), no. 5, 216-219.
4. J. L. Aurentz, T. Mach, R. Vandebril, and D. S. Watkins, *Fast and backward stable computation of roots of polynomials*, SIAM J. Matrix Anal. Appl. 36 (2015), no. 3, 942-973.
5. D. A. Bini, *Numerical computation of polynomial zeros by means of Aberth's method*, Numer. Algorithms 13 (1996), 179-200
6. E. Hansen, M. Patrick, and J. Rusnak, *Some modifications of Laguerre's method*, BIT 17 (1977), no. 4, 409-417.
7. P. Lancaster, *Lambda-matrices and vibrating systems*, International Series of Monographs on Pure and Applied Mathematics, vol. 94, Pergamon, Oxford, United Kingdom, 1966. 
8. B. Parlett, *Laguerre's method applied to the matrix eigenvalue problem*, Math. Comp. 18 (1964), no. 87, 464-485.
9. F. Tisseur, *Backward error and condition of polynomial eigenvalue problems*, Linear Algebra Appl. 309 (2000), no. 1-3, 339-361. 