# Synthetic MRI 1

Revised version.


**This code needs [Eigen](http://eigen.tuxfamily.org/) - a C++ library header files**. 
The header files location must be with the proper PATH, or PATH should be added in time of compilation.
The **optimizer** also uses Eigen, and recent versions can be found [here](https://github.com/PatWie/CppNumericalSolvers).
It uses [gsl library](https://www.gnu.org/software/gsl/) for bessel functions and [openmp](https://www.openmp.org/) for parallel processing.

Locations:
* The file to be executed for 3D ECM: ./examples/3D/example_AECM.cpp (or ./examples/3D/scheme_new_OSL_EM_29_GEM.cpp, or number 25 for single core)
    
    First go to examples/3D/
    Then compile:
        g++ example_AECM.cpp -o example_AECM -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp -DEIGEN_DONT_PARALLELIZE
    
    Then run:
        ./example_AECM ../data/ZHRTS1.nii Dummy_sd_3D.txt 0

* The file to be executed for 2D ECM: ./examples/2D/example_AECM.cpp (or ./examples/2D/scheme_new_OSL_EM_29_GEM.cpp, or number 25 - single core)
    
    First go to examples/2D
    Then compile:
        g++ example_AECM.cpp -o example_AECM -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp -DEIGEN_DONT_PARALLELIZE
    
    Then run:
        ./example_AECM ../data/new_phantom.nii Dummy_sd.txt 0


* The file to be executed for OSL-EM: 
	./examples/3D/example_OSL.cpp
	or 
	./examples/2D/example_OSL.cpp

* The data: ./data/ZHRTS1.nii (or ./data/small.nii)

	OR
	    ./data/new_phantom.nii (or ./data/small_phantom.nii)
  (see `*' for any 2D file)
* The optimizer location: ./CppNumericalSolvers


The current tree structure is as follows:
```bash
.
|-- data
|   |-- Dummy_sd.txt
|   |-- new_phantom.nii
|   |-- phantom.nii
|   |-- smallest_phantom.nii
|   |-- small.nii
|   |-- small_phantom.nii
|   `-- ZHRTS1.nii
|-- docs
|   `-- Doxyfile
|-- examples
|   |-- 2D
|   |   |-- Dummy_sd.txt
|   |   |-- example_AECM.cpp
|   |   |-- example_OSL.cpp
|   |   |-- example_VAR.cpp
|   |   `-- result
|   `-- 3D
|       |-- All_possible_AECM.cpp
|       |-- All_possible_OSL.cpp
|       |-- Dummy_sd_3D.txt
|       |-- Dummy_sd_brainweb.txt
|       |-- example_AECM.cpp
|       |-- example_OSL.cpp
|       |-- example_VAR.cpp
|       |-- example_VAR_part.cpp
|       |-- result
|       `-- Var.cpp
|-- include
|   |-- 2D
|   |   |-- functions_AECM.hpp
|   |   |-- functions_gen.hpp
|   |   |-- functions_LS_and_init_value.hpp
|   |   |-- functions_OSL.hpp
|   |   |-- functions_VAR.hpp
|   |   `-- read_files.hpp
|   `-- 3D
|       |-- functions_AECM.hpp
|       |-- functions_gen.hpp
|       |-- functions_LS_and_init_value.hpp
|       |-- functions_OSL.hpp
|       |-- functions_VAR.hpp
|       `-- read_files.hpp
|-- matlab
|-- R
|   |-- 2D_AECM.cpp
|   |-- 2D_AECM.hpp
|   |-- 2D_AECM_old.hpp
|   |-- Analyse_Hessian.R
|   `-- Analyse_Var_est.R
|-- README.md
`-- TODO

```




(* new_phantom.nii is actually transformed from phantom.nii(2D)

You can use R('oro.nifti') to read phantom.nii and then use
dim_(phantom) <- c(4, 256, 256, 1, 18, 1, 1, 1) # or equivalent
to change the dimension - as the dim are written in X, Y, Z, T/M - in this order.
It wold be directly incorporated through Read_files_2.cpp later. 

To change the dim, dim function of base R  would not work.
'RNifti' package is faster, but can't change dim of internal image.
)
