# Synthetic MRI 1

Very primitive situation in git - may contain errors.

Revised version.


**This code needs [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - a C++ library header files**. 
The header files location must be in the PATH, or should be added in time of compilation.
Th optimizer also uses Eigen, and can be found [here](https://github.com/PatWie/CppNumericalSolvers)

Locations:
* The file to be executed for ECM: ./TRY_EIGEN_CODE/scheme_new_EM_12_GEM.cpp

    To compile:
        g++ scheme_new_EM_12_GEM.cpp -o test -I /usr/include/eigen3 -O3

    To run:
        ./test ../data/new_phantom.nii ../data/Dummy_sd.txt 0

* The file to be executed for EM: ./TRY_EIGEN_CODE/scheme_new_EM_5_numerical_cholesky.cpp
* The data: ./data/ZHRTS1.nii

	OR
	    ./data/new_phantom.nii
  (see `*' for any 2D file)
* The optimizer location: ./CppNumericalSolvers
* The old optimiser location: ./optim_cpp_solver/

* Check!




(* new_phantom.nii is actually transformed from phantom.nii(2D)

You can use R('oro.nifti') to read phantom.nii and then use
dim_(phantom) <- c(4, 256, 256, 1, 18, 1, 1, 1) # or equivalent
to change the dimension - as the dim are written in X, Y, Z, T/M - in this order.
It wold be directly incorporated through Read_files_2.cpp later. 

To change the dim, dim function of base R  would not work.
'RNifti' package is faster, but can't change dim of internal image.
)
