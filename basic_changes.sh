find ./src -type f -exec sed -i -e 's/Rcpp::Rcerr/std::cerr/g' {} \;
find ./src -type f -exec sed -i -e 's/Rcpp::Rcout/std::cout/g' {} \;

