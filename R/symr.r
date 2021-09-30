#' Synthetic Magnetic Resonance (symr)
#'
#' @rdname symr
#' @param W A numeric Matrix of size n x 3, supplied as the initial value for the problem. One does not have to supply it for LS. Otherwise, if NULL, a basic estimate of W is used.
#' @param method The method to be used, possible options are.
#' @param dimen The dimension of the train MR signals (possibly read from the nifti file, the format would be: )
#' @param train A numeric matrix of size n x m of observed signals to be trained on, where m is the number of training MR images.
#' @param TE.train A numeric vector, TE values for the training set
#' @param TR.train A numeric vector, TR values for the training set
#' @param sigma.train A numeric vector, sigma_j values for the training set
#' @param train.scale A positive real number by which voxels are scaled
#' @param TE.scale A positive real number by which TE values are scaled
#' @param TR.scale A positive real number by which TR values are scaled
#' @param black.list A numeric vector of size n signifying the background voxels
#' @param maxiter.LS The maximum iteration number for the L-BFGS-B procedure for LS
#' @param maxiter The maximum iteration number for the EM algorithms for MLE/OSL/AECM
#' @param abs.diff Absolute difference criteria to stop the EM algorithm
#' @param rel.diff Relative difference criteria to stop the EM algorithm
#' @param verbose verbose outputs
#' @param verbose2 More verbose outputs
#' @return The final estimate of \code{W} (and estimate of Psi inverse and beta whenever possible) after executing the method.
#' @export
#' @examples
#' ## Basic 2D example:
#' ### Load an nifti file (using oro.nifti or Rnifti or similar package) and resizing into size n x m:
#' file_name <- system.file("extdata", "new_phantom.nii.gz", package = "symR", mustWork = TRUE)
#' phantom <- RNifti::readNifti(file_name, internal = TRUE)
#' phantom <- apply(phantom, 4, function(x) {
#'   c(x)
#' })
#' phantom[phantom == 0.0] <- 0.5 ## Pre-processing to remove the -Inf issue in likelihood.
#' n <- nrow(phantom)
#'
#' ## Other input parameters:
#' TE_values <- c(
#'   0.03, 0.06, 0.04, 0.08, 0.05, 0.10, 0.03, 0.06, 0.04,
#'   0.08, 0.05, 0.10, 0.03, 0.06, 0.04, 0.08, 0.05, 0.10
#' )
#' TR_values <- c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3)
#' sigma_values <- c(
#'   19.6926, 19.5212, 20.4545, 14.9832, 18.9208, 13.5797,
#'   21.9965, 21.4306, 25.0911, 21.1322, 22.1558, 18.9088,
#'   27.0099, 23.961, 25.0904, 22.2281, 26.9848, 22.1567
#' )
#' TE.scale <- 2.01 / min(TE_values)
#' TR.scale <- 2.01 / min(TR_values)
#' r_scale <- 10.0
#' TE_values <- TE_values * TE.scale
#' TR_values <- TR_values * TR.scale
#' phantom <- phantom / r_scale
#' sigma_values <- sigma_values / r_scale
#'
#' ## Make a mask or supply the mask:
#' mask <- 1 - {
#'   rowSums(phantom > 50) > 0
#' }
#'
#' ## Divide into train and test with 3 train images:
#' train_ind <- c(1, 7, 14)
#' test_ind <- setdiff(1:ncol(phantom), train_ind)
#' train <- phantom[, train_ind]
#' sigma.train <- sigma_values[train_ind]
#' TE.train <- TE_values[train_ind]
#' TR.train <- TR_values[train_ind]
#' test <- phantom[, test_ind]
#' sigma_test <- sigma_values[test_ind]
#' TE_test <- TE_values[test_ind]
#' TR_test <- TR_values[test_ind]
#'
#' dimen <- c(3, 256, 256, 1) ## First element correspond to dim+1
#'
#'
#' ## Get LS estimate:
#' W_LS <- symr(NULL,
#'   method = "LS", dimen, TE.train, TR.train, sigma.train, train,
#'   r_scale, TE.scale, TR.scale, mask
#' )
#'
#' ## Get AECM estimate starting from LS initial value:
#' W_AECM <- symr(W_LS,
#'   method = "AECM", dimen, TE.train, TR.train, sigma.train, train,
#'   r_scale, TE.scale, TR.scale, mask
#' )
#'
#' # Overall performaces of LS
#' mean(performance(W_LS, test, TE_test, TR_test, sigma_test, mask, 1, 1, 1))
symr <- function(W = NULL,
                 method = c("LS", "Least Square", "ML", "Maximum Likelihood", "OSL-EM", "One Step Late EM", "AECM", "EM"),
                 dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter.LS = 50L, maxiter = 50L, abs.diff = 1e-1, rel.diff = 1e-5, verbose = 0L, verbose2 = 0L) {
  if (method == "LS" | method == "Least Square") {
    W_1_init <- exp(-1 / (2.0 * TR.scale))
    W_2_init <- exp(-1 / (0.1 * TR.scale))
    .Call("_symR_Init_val_least_sq_R", PACKAGE = "symR", train, TE.train, TR.train, dimen, train.scale, TE.scale, TR.scale, maxiter.LS, W_1_init, W_2_init)
  } else {
    ##  If W is not given for other methods, take the initial simplest guess.
    if (is.null(W)) {
      n <- NROW(train)
      W <- array(1, dim = c(n, 3))
      W[, 1] <- colMeans(train)
      W[, 2] <- array(W_1_init, dim = n)
      W[, 3] <- array(W_2_init, dim = n)
      ## Numerical modifications:
      W[, 1] <- ifelse(W[, 1] > 450.0, 425.0, W[, 1])
      W[, 1] <- ifelse(W[, 1] < 0.0001, 0.0001, W[, 1])
    }
    ## So there is an estimated W now.

    ## Check 2D/3D:
    if (dimen[4] == 1) {
      ## 2D - OSL/MLE/AECM:
      if (method == "ML" | method == "Maximum Likelihood") {
        .Call("_symR_OSL_R", PACKAGE = "symR", W, dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter, 0L, abs.diff, rel.diff, verbose, verbose2)
        ## Or the AECM???? (Two step case - but more checked - Subrata)
      } else if (method == "OSL-EM" | method == "One Step Late EM") {
        .Call("_symR_OSL_R", PACKAGE = "symR", W, dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter, 1L, abs.diff, rel.diff, verbose, verbose2)
      } else if (method == "EM" | method == "AECM") {
        .Call("_symR_AECM_R", PACKAGE = "symR", W, dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter, 1L, abs.diff, rel.diff, verbose, verbose2)
      }
    } else if (dimen[4] > 1) {
      ## 3D - OSL/MLE/AECM:
      if (method == "ML" | method == "Maximum Likelihood") {
        .Call("_symR_OSL_R_3D", PACKAGE = "symR", W, dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter, 0L, abs.diff, rel.diff, verbose, verbose2)
        ## Or the AECM???? (Two step case - but more checked - Subrata)
      } else if (method == "OSL-EM" | method == "One Step Late EM") {
        .Call("_symR_OSL_R_3D", PACKAGE = "symR", W, dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter, 1L, abs.diff, rel.diff, verbose, verbose2)
      } else if (method == "EM" | method == "AECM") {
        .Call("_symR_AECM_R_3D", PACKAGE = "symR", W, dimen, TE.train, TR.train, sigma.train, train, train.scale, TE.scale, TR.scale, black.list, maxiter, 1L, abs.diff, rel.diff, verbose, verbose2)
      }
    }
  }
}
