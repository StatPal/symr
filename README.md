# Synthetic MRI 1

Very primitive situation in git - may contain errors.
Revised version.

**This code needs Eigen - a C++ library header files**. The header files location must be in the PATH, or should be added in time of compilation.


Locations:
* The file to be executed for ECM: ./TRY_EIGEN_CODE/scheme_new_EM_12_GEM.cpp
    To compile:
        g++ scheme_new_EM_12_GEM.cpp -o test -I /usr/include/eigen3 -O3
    To run:
        ./test ../data/new_phantom.nii ../data/Dummy_sd.txt 0
* The file to be executed for EM: ./TRY_EIGEN_CODE/scheme_new_EM_5_numerical_cholesky.cpp
* The data location: ./data/ZHRTS1.nii
* The optimizer location: ./CppNumericalSolvers
* The old optimiser location: ./optim_cpp_solver/

* Check!