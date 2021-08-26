# https://stackoverflow.com/questions/65872960/use-valgrind-with-r-cmd-check

#devtools::load_all()   ## This is creating problem too
#devtools::run_examples(start='symR-package.Rd')
#devtools::run_examples(fresh=TRUE)

## Somehow these give error - but the original code does not. 










library(symR)
## Basic 2D example: 
### Load an nifti file (using oro.nifti or Rnifti or similar package) and resizing into size nxm:
file_name <- system.file("extdata", "new_phantom.nii.gz", package = "symR", mustWork = TRUE)
#RNifti::readNifti("R/Researchs/MRI/Git/symR/inst/extdata/new_phantom.nii.gz")
phantom <- RNifti::readNifti(file_name, internal=TRUE)
phantom <- apply(phantom, 4, function(x){c(x)})
phantom[phantom==0.0] <- 0.5  ## Pre-processing to remove the -Inf issue in likelihood. 
n <- nrow(phantom)


## Other input parameters: 
TE_values <- c(0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 
   0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10)
TR_values <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3)
sigma_values <- c(19.6926, 19.5212, 20.4545, 14.9832, 18.9208, 13.5797, 
  21.9965, 21.4306, 25.0911, 21.1322, 22.1558, 18.9088, 
  27.0099, 23.961, 25.0904, 22.2281, 26.9848, 22.1567)
TE_scale = 2.01/min(TE_values); TR_scale = 2.01/min(TR_values); r_scale = 10.0
TE_values = TE_values*TE_scale; TR_values = TR_values*TR_scale;
phantom = phantom/r_scale; sigma_values <- sigma_values/r_scale

## Make a mask or supply the mask: 
mask <- array(1, dim=n)
for (i in 1:n) {
   mask[i] <- ifelse(any(phantom[i,]>50), 0, 1)
}

## Divide into train and test with 3 train images: 
train_ind <- c(1, 7, 14)
test_ind <- setdiff(1:ncol(phantom), train_ind)
train <- phantom[,train_ind]; sigma_train <- sigma_values[train_ind]
TE_train <- TE_values[train_ind]; TR_train <- TR_values[train_ind]
test <- phantom[,test_ind]; sigma_test <- sigma_values[test_ind]
TE_test <- TE_values[test_ind]; TR_test <- TR_values[test_ind]

dimen <- c(3, 256, 256, 1)		## First element correspond to dim+1
W1_init <- exp(-1/(2.0*TR_scale))
W2_init <- exp(-1/(0.1*TE_scale))



## Get LS estimate: 
W_init <- Init_val_least_sq_R(train, TE_train, TR_train, dimen, 
                 r_scale, TE_scale, TR_scale, W1_init, W2_init)

# Overall performaces of LS
mean(Performance_test_R(W_init, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))

## Predicted images: 
est_LS <- v_mat_R(W_init, TE_test, TR_test)
# First image with proper dimension: 
estimate_1 <- matrix(est_LS[,1], dimen[2], dimen[3])




## Get MPLE estimates using AECM:
AECM_val <- AECM_R(W_init, dimen, TE_train, TR_train, sigma_train, train, 
   r_scale, TE_scale, TR_scale, mask, 
   maxiter = 10, penalized = 1, abs_diff = 1e-1, 
   rel_diff = 1e-3)
## For possibly better results, increase maxiter/decrease rel_diff or abs_diff

# Overall performaces of AECM
mean(Performance_test_R(AECM_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))

## Predicted images: 
est_AECM <- v_mat_R(AECM_val$W, TE_test, TR_test)
# First image with proper dimension: 
estimate_1 <- matrix(est_AECM[,1], dimen[2], dimen[3])
cat("AECM done")





#### Somehow there is a problem putting it before ### CAREFULLY CHECK


## Sample row of parameters, W
W_row <- c(50, 0.01, 0.003)
W <- rbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04))	## Two sample rows

## Design parameters
TE <- c(0.01, 0.03, 0.04, 0.01)
TR <- c(0.6, 0.6, 1, 0.8)
sig <- c(1.2, 1.1, 1.4, 1.2)

## Forward transformed values: 
Bloch_eqn_R(W_row, TE, TR)
v_mat_R(W, TE, TR)
Generate_r(W, TE, TR, sig)

dee_v_ij_dee_W_ik(W_row, TE, TR, 1, 1)					## dv_i1/dW_i1
dee_2_v_ij_dee_W_ik_dee_W_ik1(W_row, TE, TR, 1, 1, 2)	## d^2v_i1/(dW_i1 dW_i2)






library(Matrix)
mean_rice_R(2,1.5)
J(5)
eigenvals_J(5)




W <- rbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04))		## Two sample rows of parameters, W
test <- rbind(c(56, 52, 57, 51), c(39, 37, 33, 34.4) )

## Test design parameters
TE <- c(0.01, 0.03, 0.04, 0.01)
TR <- c(0.6, 0.6, 1, 0.8)
sig <- c(1.2, 1.1, 1.4, 1.2)
mask <- c(0,0)


Performance_test_R(W, test, TE, TR, sig, mask, 1, 1, 1)
Performance_test_R(W, test, TE, TR, sig, mask, 3, 1, 1)
Performance_test_R(W, test, TE, TR, sig, mask, 1, 2, 1)
Performance_test_R(W, test, TE, TR, sig, mask, 3, 2, 1)











##devtools::build_vignettes()
devtools::test()
