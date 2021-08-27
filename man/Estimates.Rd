\name{Estimation-functions}
\alias{Init_val_least_sq_R}
\alias{AECM_R}
\alias{OSL_R}
\alias{AECM_R_3D}
\alias{OSL_R_3D}
\title{Set of functions for estimation of Synthetic MRI parameters}
\description{
  Set of main functions for LS/OSL/AECM estimation of Synthetic MRI parameters
}
\usage{
Init_val_least_sq_R(train, TE_train, TR_train, our_dim_1, 
                     train_scale, TE_scale, TR_scale, maxiter_LS, W_1_init, W_2_init)
AECM_R(W, our_dim_1, TE_train, TR_train, sigma_train, train, 
       train_scale, TE_scale, TR_scale, black_list, 
	   maxiter = 50, penalized = 1, abs_diff = 1e-1, 
	   rel_diff = 1e-5, verbose = 0, verbose2 = 0)
OSL_R(W, our_dim_1, TE_train, TR_train, sigma_train, train, 
       train_scale, TE_scale, TR_scale, black_list, 
	   maxiter = 50, penalized = 1, abs_diff = 1e-1, 
	   rel_diff = 1e-5, verbose = 0, verbose2 = 0)
AECM_R_3D(W, our_dim_1, TE_train, TR_train, sigma_train, train, 
          train_scale, TE_scale, TR_scale, black_list, 
	      maxiter = 50, penalized = 1, abs_diff = 1e-1, 
	      rel_diff = 1e-5, verbose = 0, verbose2 = 0)
OSL_R_3D(W, our_dim_1, TE_train, TR_train, sigma_train, train, 
         train_scale, TE_scale, TR_scale, black_list, 
	     maxiter = 50, penalized = 1, abs_diff = 1e-1, 
	     rel_diff = 1e-5, verbose = 0, verbose2 = 0)
}
\arguments{
  \item{W}{a numeric Matrix of size nx3, supplied as the initial value for the problem}
  \item{train}{a numeric Matrix, supplied as the values of the voxels of size nxm}
  \item{our_dim_1}{a numeric vector showing the dimension of the images}
  \item{TE_train}{a numeric vector, TE values for the training set}
  \item{TR_train}{a numeric vector, TR values for the training set}
  \item{sigma_train}{a numeric vector, sigma_j values for the training set}
  \item{train_scale}{a positive real number by which voxels are scaled}
  \item{TE_scale}{a positive real number by which TE values are scaled}
  \item{TR_scale}{a positive real number by which TR values are scaled}
  \item{maxiter_LS}{maximum number of iteration for LS L-BFGS-B algorithm}
  \item{black_list}{a numeric vector representing the background voxels}
  \item{maxiter}{maximum iteration number of the EM algorithm}
  \item{penalized}{whether to use the MRF penalty or not. If 0, tthe answer becomes the MLE}
  \item{abs_diff}{minimum absolute difference among the successive log likelihood tolerated}
  \item{rel_diff}{minimum relative difference among the successive log likelihood tolerated}
  \item{verbose}{More verbose outputs}
  \item{verbose2}{Check again}
  \item{W_1_init}{Redundant here - del}
  \item{W_2_init}{Redundant here - del}
}
\value{
  \code{Init_val_least_sq_R()} returns the value of Least Square estimated parameter W.

  \code{AECM_R()} returns a List of estimated parameters W, Psi and beta.
  
  \code{OSL_R()} returns a List of estimated parameters W, Psi and beta.
  
  \code{AECM_R_3D()} returns a List of estimated parameters W, Psi and beta for 3D data.
  
  \code{OSL_R_3D()} returns a List of estimated parameters W, Psi and beta for 3D data.
  
}
\details{
	These functions estimates the Synthetic MR image parameters from the train set of images which corresponds to TE_train, TR_train and sigma_train. 
}
\seealso{
  See \code{\link{symR-package}} for more details.
}
\examples{
## 
}
\author{Subrata Pal, Somak Dutta, Ranjan Maitra}

