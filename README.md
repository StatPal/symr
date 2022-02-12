# Synthetic MRI

## What is symr?

*symr* is a C++ software for Synthetic Magnetic Resonance (MR) technique which predicts images at new design parameter from few observed MR scans. The speciality of the method behind the *symr* is that it carefully uses both the physical and statistical properties of the underlying MR signal and noise. This is a theoretically sound and computationally practical matrix-free approach using multi-layered Gaussian Markov Random Field, which can predict images from as small as three MR scans, which can be used in individualized patient- and anatomy-specific contexts. We have also developed an accurate estimation of standard errors of the regional means in in the predicted images. 

*symR*, or the *R* branch(i.e., this branch with [link](https://github.com/StatPal/symr/tree/R)) consists of the R wrapper of the C++ function.  



## Installation and Usage:

### Dependencies:

Your system originally needs following softwares installed
* C++ compiler (tested with [GCC, the GNU Compiler Collection](https://gcc.gnu.org))
* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) - a header only C++ library
* [gsl library](https://www.gnu.org/software/gsl/) for Bessel functions
* [openmp](https://www.openmp.org/) for parallel processing.

We have used an **optimizer** in C++, which also uses Eigen, and recent versions can be found [here](https://github.com/PatWie/CppNumericalSolvers)

### Download:

The simplest way to install this package in R is to use [devtools](https://CRAN.R-project.org/package=devtools) package like:
```R
devtools::install_github("StatPal/symR@R")
```

Otherwise, in Linux/mac you can clone or download this branch of this repository and use the command,
```bash
R CMD build symr-R
R CMD install symR_1.0.tzr.gz
```
Without using git, you can go to the [link](https://github.com/StatPal/symr) to branch R, download it (you may directly download from this [link](https://github.com/StatPal/symr/archive/refs/heads/R.zip)) and unzip it and then issue the previous command with R installed. 
 

Additional commands to test the package (inside the directory)
```R
devtools::check()
devtools::lint()
devtools::test()
# devtools::test_coverage()
devtools::run_examples()
# devtools::check_built()
devtools::spell_check()
```

