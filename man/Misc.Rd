\name{Misc-functions}
\alias{mean_rice_R}
\alias{J}
\alias{eigenvals_J}
\title{Set of miscellaneous functions}
\description{
  Set of miscellaneous functions for estimation. 
}
\usage{
  mean_rice_R(nu, sigma)
  J(n)
  eigenvals_J(n)
}
\arguments{
  \item{nu}{Location parameter of the Rice distribution}
  \item{sigma}{Scale parameter of the Rice distribution}
  \item{n}{An integer, specifying the dimension of the J matrix}
}
\value{
  \code{mean_rice_R()} returns the value of mean of Rice distribution with location and scale parameters nu and sigma.

  \code{J()} returns the sparse matrix corresponding to Besag's Intrinsic Conditional Auto-Regressive (ICAR) matrix.
  
  \code{eigenvals_J()} returns eigenvalues corresponding to J matrix. 
}
\details{
	## 
}
\seealso{
  See Spatial interaction and the statistical analysis of lattice systems, J. Besag, 1974 for more details
}
\examples{
library(Matrix)
mean_rice_R(2,1.5)
J(5)
eigenvals_J(5)
}
\author{Subrata Pal, Somak Dutta, Ranjan Maitra}

