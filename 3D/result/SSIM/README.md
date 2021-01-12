# Estimation of the Structural Similarity Index (SSIM) for Image Quality Assessment 

Supplementary material to **Statistical estimation of the structural similarity index for image quality assessment** by Felipe Osorio, Ronny Vallejos, Wilson Barraza, Silvia Ojeda and Marcos A. Landi

Code written by: Felipe Osorio and Wilson Barraza

Correspondence author: Felipe Osorio, Email: felipe.osorios@usm.cl

### Code tested on:
- R under development (2018-02-21 r74285), running Linux Mint 18.3 (64 bits)
- R version 3.3.0, running OS X 10.13.4 (64 bits)

### Instructions: 
To create the Dynamically Loaded (DL) library, using the console prompt move to `/code/src` directory and enter:

`R CMD SHLIB -o SSIM.so *.c`

Next, copy `SSIM.so` file to the working directory (in our case to `/case_study`), and execute the commands in `case_study.R` file (optionally, you can use `reading_TIF.R` to read data from TIF images. 'tif' R package is required).

### CONTENTS:
- code/src/ : C sources for parameter estimation and computation of the structural similarity index for images (SSIM).
- code/SSIM.R : R function for computation of SSIM.
- code/SSIM_FIT.R : R function for parameter estimation of SSIM.
- data/ : reference images used in numerical experiments (available in .RDA and .PNG formats).
- README.md : this file.
