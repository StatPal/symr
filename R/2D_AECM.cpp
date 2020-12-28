#include <RcppEigen.h>
#include <RcppGSL.h>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]

//#include "../2D/read_files.hpp"
// There is some problem with char case with this in R. 
// Let's check other cases.
#include "../2D/functions_gen.hpp"
// Problem in v_mat function - when Rcpp export is on
// Problem in Gen_r_from_v_mat function - I guess due to rnd gen
// Problem goes away if Rcpp export is removed
// same for Gen_r function

#include "../2D/functions_AECM.hpp"
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


/*** R
timesTwo(42)
*/
