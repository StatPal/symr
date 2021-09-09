#' symR, a package
#'
#' The 'symR' package is an R wrapper to the C++ software symr for Synthetic Magnetic Resonance (MR) technique. Its main objective is to predict images at new design parameters(TE, TR pair) from a few observed MR scans. The specialty of our methodology is that it carefully uses both the physical and statistical properties of the underlying MR signal and noise. This method uses a theoretically sound and computationally practical, matrix-free approach modeled by a multi-layered Gaussian Markov Random Field. It can predict images from as small as three MR scans and be used in individualized patient- and anatomy-specific contexts, without the need of many a lot of data to train on. We have also developed an accurate estimation of the regional means' standard errors in the predicted images.
#' 
#' @docType package
#' @name symR
#'
#' @useDynLib symR, .registration=TRUE
NULL
#' @rawNamespace import(RcppGSL, except=c(fastLmPure, fastLm))
#' @rawNamespace import(RcppParallel, except=c(LdFlags))
#' @import Matrix
#' @importFrom Rcpp evalCpp

