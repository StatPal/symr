/* $ID: SSIM_fit.h, last updated 2018/07/27, F. Osorio */

#ifndef SSIM_FIT_H
#define SSIM_FIT_H

#include "base.h"
#include "brent.h"
#include "matrix.h"

/* terms structure (hold response and covariables) */
typedef struct TERMS_struct {
  double
    *z,         /* responses */
    *luminance, /* luminance */
    *contrast,  /* contrast */
    *structure; /* structure */
} TERMS_struct, *TERMS;

/* structure to hold model results */
typedef struct MODEL_struct {
  DIMS dm;      /* dimension data info */
  TERMS terms;  /* model terms (response and covariables) */
  int
    *pdims;     /* dimensions */
  double
    *f,         /* mean function (SSIM) */
    *x,         /* local model matrix */
    *resid,     /* residuals */
    *weights,   /* NLS weights */
    *coef,      /* coefficient estimates */
    *scale,     /* scale estimate */
    *control;   /* control settings for the estimation algorithm */
  int
    maxiter;    /* maximun number of iterations */
  double
    tolerance;  /* convergence tolerance */
} MODEL_struct, *MODEL;

/* structure to be used by the variable metric (BFGS) algorithm */
typedef struct SSIM_workdata {
  DIMS dm;
  TERMS terms;
  double phi, fnc;
  double *x, *f, *resid, *weights, *coef, *nls_der;
} SSIM_workdata, *SSIM_WKDAT;

/* structure required to evaluate the pseudo-likelihood objective function */
typedef struct PL_workdata {
  DIMS dm;
  double phi, fnc;
  double *z, *f;
} PL_workdata, *PL_WKDAT;

/* variance estimation for an heteroscedastic regression model (to be called from R) */
extern void SSIM_fit(double *, double *, double *, double *, int *, double *, double *, double *, double *);

/* for communication with gradient_test.c */
extern MODEL SSIM_init(double *, double *, double *, double *, int *, double *, double *, double *);
extern void SSIM_free(MODEL);
extern void nls_gradient(int, double *, double *, void *);
extern void update_phi(DIMS, double *, double *, double *, double);
extern double SSIM_logLik(DIMS, double *, double *, double *);
extern void mean_surface(TERMS, DIMS, double *, double *);
extern void model_frame(TERMS, DIMS, double *, double *, double *, double *, double *, int);

#endif /* SSIM_FIT_H */
