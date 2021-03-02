# Synthetic MRI

Revised version.


**This code needs [Eigen](http://eigen.tuxfamily.org/) - a C++ library header files**. 
The header files location must be with the proper PATH (or should be added in time of compilation using -I path).
The **optimizer** also uses Eigen, and recent versions can be found [here](https://github.com/PatWie/CppNumericalSolvers).
It uses [gsl library](https://www.gnu.org/software/gsl/) for bessel functions and [openmp](https://www.openmp.org/) for parallel processing.


**Instructions:**
* The data files (and the mask files) should be in Nifti format, unzipped. 
	If gziped, unzip the files: 
	```console 
	cd ./data
	gunzip new_phantom.nii.gz
	gunzip ZHRTS2.nii.gz
	gunzip ZHRTS2_class.nii.gz
	cd ..
	```


* The file to be executed (for 2D) ECM: **./examples/2D/example_AECM.cpp**
    
    First go to examples/2D
    ```console
    cd ./examples/2D/
    ``` 
    Then compile:
    ```console
    g++ example_AECM.cpp -o example_AECM -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp
    ```
    Then run:
    ```console
    ./example_AECM ../data/new_phantom.nii Dummy_sd.txt 0
    ```
    where `../data/new_phantom.nii` is the 2D data and `Dummy_sd.txt` is the file for $\sigma_j$'s generated using sigma.cpp.

* For `OSL`, everything would be similar, just the cpp file would be changed to `example_OSL.cpp` 

* For **Variance estimate** of a contrast vector(c, of size n), we have an example file with class (generated with R package mritc)
	First go to examples/2D
    ```console
    cd ./examples/2D/
    ``` 
    Then compile:
    ```console
    g++ example_VAR_part.cpp -o example_VAR_part -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm -fopenmp
    ```
    Then run:
    ```console
    ./example_VAR_part ../data/new_phantom.nii ../data/new_phantom_class.nii Dummy_sd.txt 0
    ```
    where `new_phantom_class.nii` is the file denoting class file. 



* For 3D, you have to go to `./example/3D` instead of `./example/2D` and run everything similarly with 3D data. 



* The 2D data: ./data/new_phantom.nii (see `*' for any 2D file)
  and the 3D data ./data/ZHRTS2.nii


* The optimizer location: ./CppNumericalSolvers


The current tree structure is as follows:
```bash
.
|-- data
|   |-- new_phantom.nii.gz
|   |-- phantom.nii.gz
|   |-- smallest_phantom.nii.gz
|   |-- small.nii.gz
|   `-- small_phantom.nii.gz
|-- docs
|   `-- Doxyfile
|-- examples
|   |-- 2D
|   |   |-- Dummy_sd.txt
|   |   |-- example_AECM.cpp
|   |   |-- example_OSL.cpp
|   |   |-- example_VAR.cpp
|   |   |-- example_VAR_part.cpp
|   |   |-- result
|   |   `-- sigma.cpp
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
|       `-- sigma.cpp
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
|-- R
|   |-- 2D_AECM.cpp
|   |-- 2D_AECM.hpp
|   |-- 2D_AECM_old.hpp
|   |-- Analyse_Hessian.R
|   `-- Analyse_Var_est.R
|-- README.md
`-- TODO

```



`*' 
To create the file corresponding to the $sigma_j$(rice noise parameter) for each image if they are not present,

First go to examples/2D
```console
cd ./examples/2D/
``` 
Then compile:
```console
g++ sigma.cpp -o sigma -I /usr/include/eigen3 -O3 -lgsl -lgslcblas -lm
```
Then run:
```console
./sigma ../Read_Data/new_phantom.nii Dummy_sd.txt 0
```
where `Dummy_sd.txt` is the output file containing estimated $\sigma_j$'s, i.e., the rice noise parameters. 



(* new_phantom.nii is actually transformed from phantom.nii(2D)

For any 2D data, the dimension format should be c(4, n_x, n_y, 1, m, 1, 1, 1)

You can use R('oro.nifti') to read phantom.nii and then use
dim_(phantom) <- c(4, 256, 256, 1, 18, 1, 1, 1) # or equivalent
to change the dimension - as the dim are written in X, Y, Z, T/M - in this order.
It wold be directly incorporated through Read_files_2.cpp later. 

)
