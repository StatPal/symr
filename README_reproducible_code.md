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


Step 2: 

a) Clone files for symR and run them: 
```sh
git clone -b R https://github.com/StatPal/symR --depth 1
```

b) Go to symR directory and run `image_SE_1.R`, which will run the brainweb data for noise 1%.

Options: `DEEP_LEARNING_COMPARE <- TRUE` (near line 240) instead of `DEEP_LEARNING_COMPARE <- FALSE` will add the comparison of DL and AECM and other methods. 

Step 3: Run `Variance.R`

