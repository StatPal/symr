#' Variance estimates corresponding to contrasts using Information Matrix and Delta method
#' @rdname prediction.var
#' @param  W  A numeric Matrix of size n x 3 - the estimate of primary parameters.
#' @param  Psi.inv  A numeric Matrix of size 3 x 3 - estimate of Psi inverse, corresponding to the correlation matrix between columns of the primary parameters W
#' @param  betavec  A numeric vector of size 3 - estimate of the strength of spatial factor in each directions.
#' @param  contrast  A sparse vector corresponding to the contrast of which variance is to be estimated.
#' @param  dimen The dimension of the train MR signals (possibly read from the nifti file, the format would be: )
#' @param  train A numeric matrix of size n x m of observed signals to be trained on, where m is the number of training MR images.
#' @param  TE.train A numeric vector, TE values for the training set
#' @param  TR.train A numeric vector, TR values for the training set
#' @param  sigma.train A numeric vector, sigma_j values for the training set
#' @param  TE.test A numeric vector, TE values for the test set
#' @param  TR.test A numeric vector, TR values for the test set
#' @param  sigma.test A numeric vector, sigma_j values for the test set
#' @param  train.scale A positive real number by which voxels are scaled
#' @param  TE.scale A positive real number by which TE values are scaled
#' @param  TR.scale A positive real number by which TR values are scaled
#' @param  black.list A numeric vector of size n signifying the background voxels
#' @param  test  A numeric Matrix, supplied as the values of the voxels of size n x m
#' @param  cg_maxiter  Maximum iteration number of cojugant-gradient method
#' @param  cg_tol  tolerance number of cojugant-gradient method
#' @param  penalized  This is 1 if penalized methods are used(such as "AECM", "OSL").
#' @return The vector corresponding to variance estimate
#' @export
#'
#' @examples
#' W <- rbind(c(40, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows of parameters, W
#' test <- rbind(c(56, 52, 57, 51), c(39, 37, 33, 34.4))
#'
#' ## Test design parameters
#' TE <- c(0.01, 0.03, 0.04, 0.01)
#' TR <- c(0.6, 0.6, 1, 0.8)
#' sig <- c(1.2, 1.1, 1.4, 1.2)
#' mask <- c(0, 0)
prediction.var <- function(W, Psi.inv, betavec, contrast, dimen, TE.train, TR.train, sigma.train, train, TE.test, TR.test, sigma.test, test, train.scale, TE.scale, TR.scale, black.list, cg_maxiter = 50L, cg_tol = 1e-6, penalized = 1L) {
  .Call(`_symR_Var_contrast`, W, Psi.inv, betavec, contrast, dimen, TE.train, TR.train, sigma.train, train, TE.test, TR.test, sigma.test, test, train.scale, TE.scale, TR.scale, black.list, cg_maxiter, cg_tol, penalized)
}
