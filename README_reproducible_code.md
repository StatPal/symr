To compare with the DIPSynMRI, run both Step 1 and Step 2 
with `DEEP_LEARNING_COMPARE <- TRUE` in the R file. 

To get results without Deep Learning, run only the R file 
as in Step 2 after Step 1.a, which contains a sample dataset also.

Step 3 is for Variance Estimation of ROI example. 
A contrast vector is presented as a dataset: 
phantom_1.0mm_msles1_gray_matter.mnc.gz


Step 1:
a) Clone from git
```sh
git clone https://github.com/StatPal/DeepSynMRI --depth 1
```

b) Install reticulate and related packages with the virtual environment needed for the DL package, create some necessary directories, and run the file:

```R
install.packages("reticulate")
library(reticulate)
virtualenv_create("r-reticulate")SE_example_DL.py
virtualenv_install("r-reticulate", c("torch", "torchvision", "torchaudio", "matplotlib", "scikit-image", "pandas", "joblib", "nibabel"))
dir.create("DeepSynMRI/LS-py")
dir.create("DeepSynMRI/DL_smooth")

# Run and save values
py_run_file('SE_example_DL.py')
```


Step 2: 
a) In R, to install the package run
```{R}
install.packages("devtools")
install.packages("RcppGSL") 
# sudo apt install libgsl-dev
# Or
# sudo dnf install gsl-devel
# possibly needed to run in the terminal to install gsl

devtools::install_github("StatPal/symR@R")
```

b) Clone files for symR and run them: 
```{R}
git clone -b R https://github.com/StatPal/symR --depth 1
```

c) Run `image_SE_1.R` 

Options: DEEP_LEARNING_COMPARE <- TRUE will add the comparison of DL and AECM and other methods

Step 3: Run `Variance.R`

