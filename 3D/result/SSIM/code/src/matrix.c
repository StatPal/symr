/* $ID: matrix.c, last updated 2018/07/30, F. Osorio */

#include "matrix.h"

/* basic matrix manipulations */

double
norm_sqr(double *x, int inc, int n)
{ /* sum(x * x) */
  double ans;

  ans = F77_CALL(dnrm2)(&n, x, &inc);
  return SQR(ans);
}

void
ax_plus_y(double alpha, double *x, int incx, double *y, int incy, int n)
{ /* y <- alpha * x + y */
  F77_CALL(daxpy)(&n, &alpha, x, &incx, y, &incy);
}

void
scale_vec(double alpha, double *x, int inc, int n)
{ /* x <- alpha * x */
  F77_CALL(dscal)(&n, &alpha, x, &inc);
}

void
mult_mat(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* matrix multiplication of two conformable matrices. z <- x %*% y */
  char *transx = "N", *transy = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) Calloc(xrows * ycols, double);
  F77_CALL(dgemm)(transx, transy, &xrows, &ycols, &xcols, &one, x, &ldx, y,
                  &ldy, &zero, tmp, &xrows);
  Memcpy(z, tmp, xrows * ycols);
  Free(tmp);
}

void
crossprod(double *z, double *x, int ldx, int xrows, int xcols, double *y, int ldy, int yrows, int ycols)
{ /* cross product of two given matrices. z <- t(x) %*% y */
  char *transx = "T", *transy = "N";
  double one = 1.0, zero = 0.0, *tmp = NULL;

  /* use tmp so z can be either x or y */
  tmp = (double *) Calloc(xcols * ycols, double);
  F77_CALL(dgemm)(transx, transy, &xcols, &ycols, &xrows, &one, x, &ldx, y,
                  &ldy, &zero, tmp, &xcols);
  Memcpy(z, tmp, xcols * ycols);
  Free(tmp);
}

/* DEBUG routine */

void
print_mat(char *msg, double *x, int ldx, int nrow, int ncol)
{ /* print matrix and message (used for printf debugging) */
  int i, j;

  Rprintf( "%s\n", msg);
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++)
      Rprintf( " %10.5g", x[i + j * ldx ]);
    Rprintf( "\n" );
  }
  Rprintf( "\n" );
}
