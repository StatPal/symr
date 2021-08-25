## This does the inital changes to be done
## to transform the header files C++ code to R.



## https://stackoverflow.com/questions/6758963/find-and-replace-with-sed-in-directory-and-sub-directories
find ./ -type f -exec sed -i -e 's/Rcpp::Rcerr/Rcpp::Rcerr/g' {} \;
find ./ -type f -exec sed -i -e 's/Rcpp::Rcout/Rcpp::Rcout/g' {} \;

#find ./ -type f -exec sed -i -e 's/Rcpp::stop()/Rcpp::stop()/g' {} \;
## Change back the Cppnumerical solvers README????



# Changes from C++ package:
     # -2) Put common files in include (from 2D?)
     #        functions_gen.hpp  functions_LS_and_init_value.hpp  read_files.hpp
     # -1) Give unique names to 3D and 2D rest files
     # 1) #include <Rcpp.h> after #define MAIN_HEADER
     # 2) Rcpp::Rcout instead of Rcpp::Rcout etc
     # 3) Gen_r_from_v_mat and Gen_r should be changed. Thes are still TODO. Without even srand, there is a problem with m-twister
     #    The RNG from R will not work in Cpp. Look at this and follow: https://gallery.rcpp.org/articles/random-number-generation/
     # 4)





## Ignore the read_file.hpp if possible. and read using other R functions.
