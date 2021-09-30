## This does the inital changes to be done
## to transform the header files C++ code to R.



## https://stackoverflow.com/questions/6758963/find-and-replace-with-sed-in-directory-and-sub-directories
find ./ -type f -exec sed -i -e 's/std::cerr/Rcpp::Rcerr/g' {} \;
find ./ -type f -exec sed -i -e 's/std::cout/Rcpp::Rcout/g' {} \;

#find ./ -type f -exec sed -i -e 's/Rcpp::stop()/Rcpp::stop()/g' {} \;
## Change back the Cppnumerical solvers README????



# Changes from C++ package:
     # -2) Put common files in include (from 2D?)
     #        functions_gen.hpp  functions_LS_and_init_value.hpp  read_files.hpp
     # -1) Give unique names to 3D and 2D rest files
     # 1) #include <Rcpp.h> after #define MAIN_HEADER
     # 2) Rcpp::Rcout instead of std::cout etc
     # 3) Gen_r_from_v_mat and Gen_r should be changed. Thes are still TODO. Without even srand, there is a problem with m-twister
     #    The RNG from R will not work in Cpp. Look at this and follow: https://gallery.rcpp.org/articles/random-number-generation/
     # 4)


# Remove the whole function 'Preprocess_data' from the src/include/functions_LS_and_init_value.hpp

# Then exclude the read_files.hpp ...
## That contains most number of exit(EXIT_FAILURE)




grep -iR 'exit' src/include/



## Ignore the read_file.hpp if possible. and read using other R functions



Wed Sep 29 07:37:13 PM CDT 2021

Change the Header protector of the 3D var file.


Copy the sigma.cpp file from examples/3D to src/include 
-- convert it into an hpp file with no main and with protector



Add a //[[Rcpp::export]] in that file and just include that header file in AECM_wrapper.cpp   -- This does not work. 
Just do the usual way and replace a EXIT_FAILURE


Close the Debug options


