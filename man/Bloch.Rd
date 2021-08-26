\name{Bloch equation functions}
\alias{Bloch_eqn_R}
\alias{v_mat_R}
\alias{Generate_r}
\alias{dee_v_ij_dee_W_ik}
\alias{dee_2_v_ij_dee_W_ik_dee_W_ik1}
\title{Set of functions for estimation of Synthetic MRI parameters}
\description{
  Functions corresponding to Bloch equation and generation of synthetic data from specified parameters and their derivatives. 
}
\usage{
  Bloch_eqn_R(W_row, TE, TR)
  v_mat_R(W, TE, TR)
  Generate_r(W, TE, TR, sigma)
  dee_v_ij_dee_W_ik(W_row, TE, TR, j, k)
  dee_2_v_ij_dee_W_ik_dee_W_ik1(W_row, TE, TR, j, k, k1)
}
\arguments{
  \item{W}{a numeric Matrix of size nx3 corresponding to the transformed parameters, \eqn{\rho, \mathrm{T}_1, \mathrm{T}_2}}
  \item{W_row}{a numeric Vector of size 3, any row of W}
  \item{TE}{a numeric vector, TE values for the training set}
  \item{TR}{a numeric vector, TR values for the training set}
  \item{sigma}{a numeric vector, \eqn{\sigma_j} values for the training set for all \eqn{j}}
  \item{j}{index of image number or corresponding settings number}
  \item{k, k1}{co-ordinate of the W(i.e., \eqn{W_1} or \eqn{W_2} or \eqn{W_3} ) with respect to which the derivative is taken}
}
\value{
  \code{Bloch_eqn_R()} returns the vector corresponding to the MR signal values Bloch equation predicts(i.e., \eqn{\hat{\nu}}).
  \code{v_mat_R()} returns a Matrix corresponding to the MR signal values Bloch equation predicts for a whole image.
  \code{Generate_r()} returns a Matrix corresponding to the MR signal values Bloch equation predicts for a whole image with Ricean noise. 
  \code{dee_v_ij_dee_W_ik()} retruns the derivative of \eqn{v_{ij}} with respect to(w.r.t.) \eqn{W_{ik}}.
  \code{dee_2_v_ij_dee_W_ik_dee_W_ik1()} retruns the double derivative of \eqn{v_{ij}} w.r.t. \eqn{W_{ik}} and \eqn{W_{ik_1}}.
}
\details{
	We have used the approximate Bloch equation \eqn{\nu=\rho \exp\left(-\frac{TE_j}{T_{2i}}\right)\left\{1 - \exp\left(-\frac{TR_j}{T_{1i}}\right)\right\}}
	We have used the transformed set of variables for the \eqn{T_1, T_2}: 
	First 3 functions functions gives noisy or noise-free MR signals obtained from a specific set of W parameters using Bloch equations. 
	More specifically, the first function returns the estimated signal \eqn{\nu} for a singe pixel/voxel, according to the Bloch equation. 
	The second one returns the \eqn{\nu} for all pixels/voxels at once. 
	The third one, adds i.i.d. Ricean errors with each voxels and returns the noise-contaminated MR-signal. 
	The intensity of the parameters are determined by \eqn{\sigma_j}'s. 

	Last two functions derives the derivatives and double derivatives of the \eqn{\nu_{ij}} w.r.t. \eqn{W_{ik}} and \eqn{W_{ik_1}}, needed for miscellaneous computations.
}
\seealso{
  See \code{\link{symR-package}} for more details.
}
\examples{
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
\dontrun{Generate_r(W, TE, TR, sig)}

dee_v_ij_dee_W_ik(W_row, TE, TR, 1, 1)					## dv_i1/dW_i1
dee_2_v_ij_dee_W_ik_dee_W_ik1(W_row, TE, TR, 1, 1, 2)	## d^2v_i1/(dW_i1 dW_i2)
}
\author{Subrata Pal, Somak Dutta, Ranjan Maitra}
