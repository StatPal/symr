\name{Prediction functions}
\alias{Performance_test_R}
\title{Function to measure performances of predicted cases}
\description{
  Function corresponding to measure Performances of the predictions
}
\usage{
  Performance_test_R(W, test, TE_test, TR_test, sigma_test, black_list,
					 v_type, measure_type, scale, verbose)
}
\arguments{
  \item{W}{a numeric Matrix of size nx3, supplied as the initial value for the problem}
  \item{test}{a numeric Matrix, supplied as the values of the voxels of size nxm}
  \item{TE_test}{a numeric vector, TE values for the testing set}
  \item{TR_test}{a numeric vector, TR values for the testing set}
  \item{sigma_test}{a numeric vector, sigma_j values for the testing set}
  \item{black_list}{a numeric vector representing the background voxels}
  \item{v_type}{1 corresponds to compared with nu, and 3 corresponds to compared with the Rice mean}
  \item{measure_type}{1-abs deviation, 2-squared deviation from the mean}
  \item{scale}{Scaled measure if 1}
  \item{verbose}{More verboseity if 1}
}
\value{
  \code{Performance_test_R()} returns the vector corresponding to performance measures.
}
\seealso{
  See \code{\link{symR-package}} for more details.
}
\examples{
W <- cbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04))		## Two sample rows of parameters, W
test <- cbind(c(56, 52, 57, 51), c(39, 37, 33, 34.4) )

## Test design parameters
TE <- c(0.01, 0.03, 0.04, 0.01)
TR <- c(0.6, 0.6, 1, 0.8)
sig <- c(1.2, 1.1, 1.4, 1.2)
mask <- c(1,1)

Performance_test_R(W, test, TE, TR, sig, mask, 1, 1, 1)

}
\author{Subrata Pal, Somak Dutta, Ranjan Maitra}
