/* $ID: gradient_test.c, last updated 2018/07/27, F. Osorio */

#include "gradient_test.h"

/* static functions.. */
static void gradient_multiplicative(MODEL, double *, double *, double *);
/* ..end declarations */

void
gradient_test(double *z, double *luminance, double *contrast, double *structure,
  int *pdims, double *coef, double *scale, double *control, double *score, double *stat,
  double *logLik)
{
  MODEL model;

  model = SSIM_init(z, luminance, contrast, structure, pdims, coef, scale, control);
  gradient_multiplicative(model, score, stat, logLik);
  SSIM_free(model);
}

static void
gradient_multiplicative(MODEL model, double *score, double *stat, double *logLik)
{
  SSIM_WKDAT workdata;
  DIMS dm = model->dm;
  TERMS terms = model->terms;
  double accum = 0.0, tol = R_pow(model->tolerance, 2./3.), *coef;

  /* initialization */
  coef     = (double *) Calloc(dm->p, double);
  workdata = (SSIM_WKDAT) Calloc(1, SSIM_workdata);
  for (int i = 0; i < dm->p; i++) {
    coef[i] = (model->coef)[i];
    (model->coef)[i] = 1.0;
  }

  mean_surface(terms, dm, model->coef, model->f);
  update_phi(dm, model->scale, terms->z, model->f, tol);

  /* construct nls workdata */
  workdata->dm      = dm;
  workdata->terms   = terms;
  workdata->f       = model->f;
  workdata->x       = model->x;
  workdata->resid   = model->resid;
  workdata->weights = model->weights;
  workdata->coef    = model->coef;
  workdata->nls_der = score;
  Memcpy(&(workdata->phi), model->scale, 1); /* FIXME */

  nls_gradient(dm->p, model->coef, score, workdata);

  *logLik = SSIM_logLik(dm, model->resid, model->weights, model->scale);

  /* computing gradient statistic */
  for (int i = 0; i < dm->p; i++) {
    score[i] *= -1.0;
    accum += score[i] * (coef[i] - 1.0);
  }
  *stat = accum;

  Free(coef); Free(workdata);
}
