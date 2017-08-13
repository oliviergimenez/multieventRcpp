# Multievent capture-recapture models with Rcpp

Fit a simple multievent capture-recapture model to data using Rcpp and RcppArmadillo. 

Download the multievent.cpp and titi2.txt files and type Rcpp::sourceCpp('multievent.cpp') in the R console. This code implements the matrix product by hand. The Rcpp code is 10 times faster than basic R.

Following the advice of both Dirk Eddelbuettel and Romain Francois, I have switched to RcppArmadillo to rely on 'the matrix code by professionals and decades of testing - LAPACK and BLAS'. Now the RcppArmadillo code is 60 times faster than basic R!

The results are in the multieventrcpp.pdf file.
