/* $ID: stats.h, last updated 2018/07/23, F. Osorio */

#ifndef STATS_H
#define STATS_H

#include "base.h"

/* routines for the computations of sample statistics */
void mean_and_var(double *, int, double *, double *);
void online_covariance(double *, double *, int, double *, double *, double *, double *, double *);

#endif /* STATS_H */
