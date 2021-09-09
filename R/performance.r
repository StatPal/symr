#' Performance measure computation
#' @rdname performance
#' @param  W  A numeric Matrix of size n x 3, supplied as the initial value for the problem
#' @param  test  A numeric Matrix, supplied as the values of the voxels of size n x m
#' @param  TE  A numeric vector, TE values for the testing set
#' @param  TR  A numeric vector, TR values for the testing set
#' @param  sigma  A numeric vector, sigma_j values for the testing set
#' @param  black_list  A numeric vector representing the background voxels
#' @param  v_type  1 corresponds to compared with nu, and 3 corresponds to compared with the Rice mean
#' @param  measure_type  1-abs deviation, 2-squared deviation from the mean
#' @param  scale  Scaled measure if 1
#' @param  verbose  More verboseity if 1
#' @return The vector corresponding to performance measures
#' @export
#' 
#' @examples 
#' W <- rbind(c(40, 0.01, 0.003), c(36, 0.02, 0.04))		## Two sample rows of parameters, W
#' test <- rbind(c(56, 52, 57, 51), c(39, 37, 33, 34.4) )
#'
#' ## Test design parameters
#' TE <- c(0.01, 0.03, 0.04, 0.01)
#' TR <- c(0.6, 0.6, 1, 0.8)
#' sig <- c(1.2, 1.1, 1.4, 1.2)
#' mask <- c(0,0)
#'
#' performance(W, test, TE, TR, sig, mask, 1, 1, 1)
#' performance(W, test, TE, TR, sig, mask, 3, 1, 1)
#' performance(W, test, TE, TR, sig, mask, 1, 2, 1)
#' performance(W, test, TE, TR, sig, mask, 3, 2, 1)
performance <- function(W, test, TE, TR, sigma, black.list, v_type = 1L, measure_type = 1L, scale = 1L, verbose = 0L) {
    .Call('_symR_Performance_test_R', PACKAGE = 'symR', W, test, TE, TR, sigma, black.list, v_type, measure_type, scale, verbose)
}

