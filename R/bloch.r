#' Forward Bloch Equation (for a single voxel):
#' @rdname bloch
#' @param w A numeric Vector of size 3
#' @param TE A numeric vector, TE values for the training set
#' @param TR A numeric vector, TR values for the training set
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
#' @param TE A numeric vector, TE values for the training set
#' @param TR A numeric vector, TR values for the training set
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
