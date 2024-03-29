# Synthetic MRI

## What is symr?

*symr* is C++ software for Synthetic Magnetic Resonance (MR) technique which predicts images at new design parameter from few observed MR scans. The speciality of the method behind the *symr* is that it carefully uses both the physical and statistical properties of the underlying MR ssignal and noise. This is a theoretically sound and computationally practical matrix-free approach using multi-layered Gausssain Markov Random Field, which can predict images from as small as three MR scans, which can be used in individualized patient- and anatomy-specific contexts. We have also developed an accurate estimation of standard errors of the regional means in in the predicted images. 



## Installation and Usage:

### Dependencies:

Make sure your system have the following softwares installed
* C++ compiler (tested with [GCC, the GNU Compiler Collection](https://gcc.gnu.org))
* [Eigen](http://eigen.tuxfamily.org/) - a header only C++ library
* [gsl library](https://www.gnu.org/software/gsl/) for bessel functions
* [openmp](https://www.openmp.org/) for parallel processing.

We have used an **optimizer** in C++, which also uses Eigen, and recent versions can be found [here](https://github.com/PatWie/CppNumericalSolvers)

### Download:

As this library is header only, you have to first clone/download the repository, go to the directory and then compile and run the corresponding files.
If you have git in your system, go to a working folder and run:
```console
git clone --depth=1 https://github.com/StatPal/symr.git
cd symr
```
Without using git, you can go to the [link](https://github.com/StatPal/symr), download it and unzip it and then go to the `symr` directory/folder. 



### Usage and Instructions:

* The data files (and the mask files) should be in Nifti format, unzipped.
	If gziped, unzip the files:
	```console
	cd ./data
	gunzip new_phantom.nii.gz
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
	./example_AECM ../../data/new_phantom.nii Dummy_sd.txt 0
	```
	where `../data/new_phantom.nii` is the **2D data** and `Dummy_sd.txt` is the file for &sigma;<sub>j</sub>'s (rice noise parameter) for each image generated using sigma.cpp. (See this [subsection](#generation-of-sigmas) for details).


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
	where `../data/new_phantom_class.nii` is the file denoting class file. 



* For 3D, you have to go to `./example/3D` instead of `./example/2D` and run everything similarly with 3D data. 



* The location of 2D data: ./data/new_phantom.nii (see `**' for any 2D file)
  and the 3D data ./data/ZHRTS1.nii


* The 3rd-party optimizer location: ./CppNumericalSolvers


The current tree structure is as follows:
```bash
.
|-- data
|   `-- new_phantom.nii.gz
|-- docs
|   `-- Doxyfile
|-- examples
|   |-- 2D
|   |   |-- example_AECM.cpp
|   |   |-- example_OSL.cpp
|   |   |-- example_VAR.cpp
|   |   |-- example_VAR_part.cpp
|   |   |-- result
|   |   `-- sigma.cpp
|   `-- 3D
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
`-- README.md

```


#### Generation of sigmas
To create the file corresponding to the &sigma;<sub>j</sub>(rice noise parameter) for each image if they are not present,

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
./sigma ../../data/new_phantom.nii Dummy_sd.txt 0
```
where `Dummy_sd.txt` is the output file containing estimated &sigma;<sub>j</sub>'s, i.e., the rice noise parameters. 




(`**' new_phantom.nii is actually transformed from phantom.nii(2D)

For any 2D data, the dimension format should be c(4, n_x, n_y, 1, m, 1, 1, 1)

You can use R('oro.nifti') to read phantom.nii and then use
dim_(phantom) <- c(4, 256, 256, 1, 18, 1, 1, 1) # or equivalent
to change the dimension - as the dim are written in X, Y, Z, T/M - in this order.
It wold be directly incorporated through Read_files_2.cpp later. 
)
