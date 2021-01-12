/* $ID: stats.c, last updated 2018/07/23, F. Osorio */

#include "stats.h"

void
mean_and_var(double *x, int nobs, double *mean, double *var)
{ /* computes the sample mean and variance using an online algorithm */
  int i, n = 0;
  double accum = 0.0, diff;

  *mean = 0.0;
  for (i = 0; i < nobs; i++) {
    n++;
    diff = x[i] - *mean;
    *mean += diff / n;
    accum += diff * (x[i] - *mean);
  }
  *var = accum / n;
}

void
online_covariance(double *x, double *y, int nobs, double *xbar, double *ybar,
  double *xvar, double *yvar, double *cov)
{ /* computes the sample covariance using an online algorithm */
  int i, n = 0;
  double accum = 0.0, acc_x = 0.0, acc_y = 0.0, dx, dy;

  *xbar = *ybar = 0.0;
  for (i = 0; i < nobs; i++) {
    n++;
    dx = x[i] - *xbar;
    dy = y[i] - *ybar;
    *xbar += dx / n;
    *ybar += dy / n;
    acc_x += dx * (x[i] - *xbar);
    acc_y += dy * (y[i] - *ybar);
    accum += (n - 1) * (dx / n) * (dy / n) - accum / n;
  }
  *xvar = acc_x / n;
  *yvar = acc_y / n;
  *cov  = accum;
}
