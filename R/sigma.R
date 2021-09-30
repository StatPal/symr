#' Estimating rice noise parameter from a MR image using mixture of rice distribution.
#' @rdname estimate.sigma.j
#' @param r.col  An MR image in a vector format
#' @param min_grp  Minimum number of possible mixture groups for the Rice-mixture distribution
#' @param max_grp  Maximum number of possible mixture groups for the Rice-mixture distribution
#' @param init_iter  Maximum iteration number for initial value estimation
#' @param EM_iter  Maximum iteration number for EM estimation
#' @param eps  Relative error stopping criteria for the EM algorithm
#' @return The estimate of rice noise parameter for the image \code{r.col}
#' @export
#'
#' @examples
#' \donttest{
#' ## Sample row of image,
#' file_name <- system.file("extdata", "new_phantom.nii.gz", package = "symR", mustWork = TRUE)
#' phantom <- RNifti::readNifti(file_name, internal = TRUE)
#' phantom <- apply(phantom, 4, function(x) {
#'   c(x)
#' })
#' phantom[phantom == 0.0] <- 0.5 ## Pre-processing to remove the -Inf issue in likelihood.
#' estimate.sigma.j(phantom[, 1])
#' }
estimate.sigma.j <- function(r.col, min_grp = 5L, max_grp = 15L, init_iter = 4L, EM_iter = 15L, eps = 0.0001) {
  .Call(`_symR_est_sigma_j`, r.col, min_grp, max_grp, init_iter, EM_iter, eps)
}
