/* $ID: similarity.c, last updated 2018/07/23, F. Osorio */

#include "similarity.h"

void
SSIM(double *x, double *y, double *pars, double *eps, double *stats, double *comp)
{ /* compute structural similarity index for two images */
  long int nobs;
  double xbar, ybar, xvar, yvar, cov;
  double c1, c2, c3;
  double alpha, beta, gamma;
  double luminance, contrast, structure, ssim;

  /* get constants and parameters */
  nobs = (int) eps[0];
  c1   = eps[1];
  c2   = eps[2];
  c3   = eps[3];
  alpha = pars[0];
  beta  = pars[1];
  gamma = pars[2];

  online_covariance(x, y, nobs, &xbar, &ybar, &xvar, &yvar, &cov);

  luminance = (2. * xbar * ybar + c1) / (SQR(xbar) + SQR(ybar) + c1);
  contrast  = (2. * sqrt(xvar) * sqrt(yvar) + c2) / (xvar + yvar + c2);
  structure = (cov + c3) / (sqrt(xvar) * sqrt(yvar) + c3);
  ssim = alpha * log(luminance) + beta * log(contrast) + gamma * log(structure);

  /* save results */
  stats[0] = xbar;
  stats[1] = ybar;
  stats[2] = xvar;
  stats[3] = yvar;
  stats[4] = cov;
  comp[0]  = exp(ssim);
  comp[1]  = luminance;
  comp[2]  = contrast;
  comp[3]  = structure;
}
