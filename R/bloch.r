#' Forward Bloch Equation (for a single voxel):
#' @rdname bloch
#' @param w A numeric Vector of size 3
#' @param TE A numeric vector, TE values
#' @param TR A numeric vector, TR values
#' @return The vector corresponding to the MR signal values Bloch equation predicts(i.e., \eqn{\hat{\nu}})
#' @export
#'
#' @examples
#' ## Sample row of parameters, w
#' w <- c(50, 0.01, 0.003)
#' ## Design parameters
#' TE <- c(0.01, 0.03, 0.04, 0.01)
#' TR <- c(0.6, 0.6, 1, 0.8)
#' ## Forward transformed values:
#' bloch(w, TE, TR)
bloch <- function(w, TE, TR) {
  .Call("_symR_Bloch_eqn_R", PACKAGE = "symR", w, TE, TR)
}


#' Forward Bloch Equation (for whole MR image):
#' @rdname bloch.image
#' @param W A numeric matrix of size n x 3
#' @param TE A numeric vector, TE values
#' @param TR A numeric vector, TR values
#' @return The matrix corresponding to the MR signal values Bloch equation predicts(i.e., \eqn{\hat{\nu}})
#' @export
#'
#' @examples
#' ## Sample row of parameters, W
#' W <- rbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows
#' ## Design parameters
#' TE <- c(0.01, 0.03, 0.04, 0.01)
#' TR <- c(0.6, 0.6, 1, 0.8)
#' ## Forward transformed values:
#' bloch.image(W, TE, TR)
bloch.image <- function(W, TE, TR) {
  t(apply(W, 1, bloch, TE = TE, TR = TR))
}


#' Forward Bloch Equation (for whole MR image):
#' @rdname Generate.noisy.image
#' @param W A numeric matrix of size n x 3
#' @param TE A numeric vector, TE values
#' @param TR A numeric vector, TR values
#' @param sigma A numeric vector, Rice noise parameter values
#' @return The matrix corresponding to the Rice noise contaminated MR signal values Bloch equation predicts
#' @export
#'
#' @examples
#' ## Sample row of parameters, W
#' W <- rbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows
#' ## Design parameters
#' TE <- c(0.01, 0.03, 0.04, 0.01)
#' TR <- c(0.6, 0.6, 1, 0.8)
#' sig <- c(1.2, 2.1, 1.8)
#' ## Forward transformed values:
#' Generate.noisy.image(W, TE, TR, sigma)
Generate.noisy.image <- function(W, TE, TR, sigma) {
    .Call(`_symR_Generate_r`, W, TE, TR, sigma)
}

