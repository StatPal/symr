\name{symR-package}
\alias{symR-package}
\alias{symR}
\docType{package}
\title{
\packageTitle{symR}
}
\description{
\packageDescription{symR}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{symR}
\packageIndices{symR}
     symR is a R wrapper to the C++ software symr for Synthetic Magnetic Resonance (MR) technique which predicts images at new design parameter from few observed MR scans. The speciality of the method behind the symR is that it carefully uses both the physical and statistical properties of the underlying MR signal and noise. This is a theoretically sound and computationally practical matrix-free approach using multi-layered Gausssain Markov Random Field, which can predict images from as small as three MR scans, which can be used in individualized patient- and anatomy-specific contexts. We have also developed an accurate estimation of standard errors of the regional means in in the predicted images.
}
\author{
\packageAuthor{symR}

Maintainer: \packageMaintainer{symR}
}
\references{
	For background information, look at, \href{https://arxiv.org/abs/2103.01532}{Arxiv version of the paper} and \href{https://github.com/StatPal/symr}{C++ package}
}
\concept{ Synthetic-MRI, MRI, spatial, multivariate }
\examples{
## 2D example: 
### Load an nifti file (using oro.nifti or Rnifti or similar package) and resizing into size nxm:
file_name <- system.file("extdata", "new_phantom.nii.gz", package = "symR", mustWork = TRUE)
#oro.nifti::readNIfTI("R/Researchs/MRI/Git/symR/inst/extdata/new_phantom.nii.gz")
phantom <- oro.nifti::readNIfTI(file_name)
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

cat("\nOverall performaces of LS: \n")
mean(Performance_test_R(W_init, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))
mean(Performance_test_R(W_init, test, TE_test, TR_test, sigma_test, mask, 3, 1, 1))
mean(Performance_test_R(W_init, test, TE_test, TR_test, sigma_test, mask, 1, 2, 1))
mean(Performance_test_R(W_init, test, TE_test, TR_test, sigma_test, mask, 3, 2, 1))

## Plot the predicted images: 
est_LS <- v_mat_R(W_init, TE_test, TR_test)
for(i in 1:length(TE_test)){
  tmp <- matrix(est_LS[,i], dimen[2], dimen[3])
  image(tmp)
}


## Get MLE estimates: 
MLE_val <- AECM_R(W_init, dimen, TE_train, TR_train, sigma_train, train, 
   r_scale, TE_scale, TR_scale, mask, 
   maxiter = 50, penalized = 0, abs_diff = 1e-1, 
   rel_diff = 1e-5, verbose = 0, verbose2 = 0)

cat("\nOverall performaces of MLE: \n")
mean(Performance_test_R(MLE_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))
mean(Performance_test_R(MLE_val$W, test, TE_test, TR_test, sigma_test, mask, 3, 1, 1))
mean(Performance_test_R(MLE_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 2, 1))
mean(Performance_test_R(MLE_val$W, test, TE_test, TR_test, sigma_test, mask, 3, 2, 1))

## Plot the predicted images: 
est_MLE <- v_mat_R(MLE_val$W, TE_test, TR_test)
for(i in 1:length(TE_test)){
  tmp <- matrix(est_MLE[,i], dimen[2], dimen[3])
  image(tmp)
}



# Get MPLE estimates using OSL: 
OSL_val <- OSL_R(W_init, dimen, TE_train, TR_train, sigma_train, train, 
   r_scale, TE_scale, TR_scale, mask, 
  maxiter = 50, penalized = 1, abs_diff = 1e-1, 
  rel_diff = 1e-5, verbose = 0, verbose2 = 0)

cat("\nOverall performaces of OSL: \n")
mean(Performance_test_R(OSL_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))
mean(Performance_test_R(OSL_val$W, test, TE_test, TR_test, sigma_test, mask, 3, 1, 1))
mean(Performance_test_R(OSL_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 2, 1))
mean(Performance_test_R(OSL_val$W, test, TE_test, TR_test, sigma_test, mask, 3, 2, 1))

## Plot the predicted images: 
est_OSL <- v_mat_R(OSL_val$W, TE_test, TR_test)
for(i in 1:length(TE_test)){
  tmp <- matrix(est_OSL[,i], dimen[2], dimen[3])
  image(tmp)
}



## Get MPLE estimates using AECM:
AECM_val <- AECM_R(W_init, dimen, TE_train, TR_train, sigma_train, train, 
   r_scale, TE_scale, TR_scale, mask, 
   maxiter = 50, penalized = 1, abs_diff = 1e-1, 
   rel_diff = 1e-5, verbose = 0, verbose2 = 0)
# head(AECM_val$W); AECM_val$Psi_inv; AECM_val$beta

cat("\nOverall performaces of AECM: \n")
mean(Performance_test_R(AECM_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))
mean(Performance_test_R(AECM_val$W, test, TE_test, TR_test, sigma_test, mask, 3, 1, 1))
mean(Performance_test_R(AECM_val$W, test, TE_test, TR_test, sigma_test, mask, 1, 2, 1))
mean(Performance_test_R(AECM_val$W, test, TE_test, TR_test, sigma_test, mask, 3, 2, 1))

## Plot the predicted images: 
est_AECM <- v_mat_R(AECM_val$W, TE_test, TR_test)
for(i in 1:length(TE_test)){
  tmp <- matrix(est_AECM[,i], dimen[2], dimen[3])
  image(tmp)
}



## 3D example: 
}
