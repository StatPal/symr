/* $ID: base.h, last updated 2018/07/30, F. Osorio */

#ifndef BASE_H
#define BASE_H

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>

/* some definitions */
#define DNULLP      (double *) 0
#define MAX(a,b)    (((a)>(b)) ? (a) : (b))
#define MIN(a,b)    (((a)<(b)) ? (a) : (b))
#define EQUAL(a,b)  (((a)!=(b)) ? (0) : (1))
#define SQR(x)      R_pow_di(x, 2)
#define ABSTOL      1.0e-5
#define REPORT      1
#define GOLDEN      0.3819660112501051
#define repeat      for(;;)

/* dims structure */
typedef struct DIMS_struct {
  int
    n,  /* number of observations */
    p;  /* number of parameters */
} DIMS_struct, *DIMS;

#endif /* BASE_H */
