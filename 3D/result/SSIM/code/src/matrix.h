/* $ID: matrix.h, last updated 2018/07/30, F. Osorio */

#ifndef MATRIX_H
#define MATRIX_H

#include "base.h"

/* basic matrix manipulations and BLAS wrappers */
extern double norm_sqr(double *, int, int);
extern void ax_plus_y(double, double *, int, double *, int, int);
extern void scale_vec(double, double *, int, int);
extern void mult_mat(double *, double *, int, int, int, double *, int, int, int);
extern void crossprod(double *, double *, int, int, int, double *, int, int, int);

/* DEBUG routine */
extern void print_mat(char *, double *, int, int, int);

#endif /* MATRIX_H */
