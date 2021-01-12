/* $ID: SSIM_fit.c, last updated 2018/07/27, F. Osorio */

#include "SSIM_fit.h"

/* static functions.. */
static DIMS dims(int *);
static void dims_free(DIMS);
static TERMS terms_init(double *, double *, double *, double *);
static void terms_free(TERMS);
static int SSIM_iterate(MODEL);
static void update_coef(MODEL);
static double nls_objective(int, double *, void *);
static void objective_SSIM(int, double *, double *, double *, int, SSIM_WKDAT);
static double PL_objective(double, void *);
/* ..end declarations */

void
SSIM_fit(double *z, double *luminance, double *contrast, double *structure, int *pdims,
  double *coef, double *scale, double *control, double *logLik)
{ /* SSIM parameter estimation using a heteroscedastic regression model */
  MODEL model;

  model = SSIM_init(z, luminance, contrast, structure, pdims, coef, scale, control);
  control[2] = (double) SSIM_iterate(model);
  *logLik = SSIM_logLik(model->dm, model->resid, model->weights, model->scale);
  SSIM_free(model);
}

static DIMS
dims(int *pdims)
{ /* constructor for a dims object */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->n = (int) pdims[0];
  ans->p = (int) pdims[1];
  return ans;
}

static void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}

static TERMS
terms_init(double *z, double *luminance, double *contrast, double *structure)
{ /* constructor for a terms object */
  TERMS ans;

  ans = (TERMS) Calloc(1, TERMS_struct);
  ans->z = z;
  ans->luminance = luminance;
  ans->contrast  = contrast;
  ans->structure = structure;
  return ans;
}

static void
terms_free(TERMS this)
{ /* destructor for a terms object */
  Free(this);
}

MODEL
SSIM_init(double *z, double *luminance, double *contrast, double *structure, int *pdims,
  double *coef, double *scale, double *control)
{ /* constructor for a heteroscedastic regression model object */
  MODEL model;

  model = (MODEL) Calloc(1, MODEL_struct);
  model->dm = dims(pdims);
  model->terms = terms_init(z, luminance, contrast, structure);
  model->f = (double *) Calloc(model->dm->n, double);
  model->x = (double *) Calloc(model->dm->n * model->dm->p, double);
  model->resid   = (double *) Calloc(model->dm->n, double);
  model->weights = (double *) Calloc(model->dm->n, double);
  model->coef  = coef;
  model->scale = scale;
  model->control = control;
  model->maxiter = (int) control[0];
  model->tolerance = control[1];
  return model;
}

void
SSIM_free(MODEL this)
{ /* destructor for a model object */
  dims_free(this->dm);
  terms_free(this->terms);
  Free(this->f);
  Free(this->x);
  Free(this->resid);
  Free(this->weights);
  Free(this);
}

static int
SSIM_iterate(MODEL model)
{ /* iteratively weighted NLS algorithm */
  int iter = 0, job;
  double conv, logLik, newlogLik, tol = R_pow(model->tolerance, 2./3.);

  /* initialization */
  job = 0;
  model_frame(model->terms, model->dm, model->coef, model->scale, model->x,
              model->resid, model->weights, job);
  logLik = SSIM_logLik(model->dm, model->resid, model->weights, model->scale);

  /* main loop */
  repeat {
    /* PL-step */
    mean_surface(model->terms, model->dm, model->coef, model->f);
    update_phi(model->dm, model->scale, model->terms->z, model->f, tol);

    /* NLS-step */
    job = 1;
    model_frame(model->terms, model->dm, model->coef, model->scale, model->x,
                model->resid, model->weights, job);
    update_coef(model);
    newlogLik = SSIM_logLik(model->dm, model->resid, model->weights, model->scale);

    iter++;

    /* eval convergence */
    conv = fabs((newlogLik - logLik) / (newlogLik + ABSTOL));
    if (conv < tol)
      break; /* successful completion */
    if (iter >= model->maxiter)
      break; /* maximum number of iterations exceeded */
    logLik = newlogLik;
  }
  return iter;
}

static void
update_coef(MODEL model)
{
  SSIM_WKDAT workdata;
  DIMS dm = model->dm;
  TERMS terms = model->terms;
  int fail = 0, fncount, grcount, trace = 0, *mask;
  double *coef, *gradient, nls_fnc = 0.0;
  char *vmax;

  workdata = (SSIM_WKDAT) Calloc(1, SSIM_workdata);
  mask     = (int *) Calloc(dm->p, int);
  coef     = (double *) Calloc(dm->p, double);
  gradient = (double *) Calloc(dm->p, double);

  for (int i = 0; i < dm->p; i++) {
    coef[i] = (model->coef)[i];
    mask[i] = 1;
  }

  /* construct nls workdata */
  workdata->dm      = dm;
  workdata->terms   = terms;
  workdata->f       = model->f;
  workdata->x       = model->x;
  workdata->resid   = model->resid;
  workdata->weights = model->weights;
  workdata->coef    = coef;
  workdata->nls_der = gradient;
  Memcpy(&(workdata->phi), model->scale, 1); /* FIXME */

  /* Calling variable metric (BFGS) optimizer */
  vmax = vmaxget();
  vmmin(dm->p, coef, &nls_fnc, nls_objective, nls_gradient, model->maxiter, trace,
        mask, (double) ABSTOL, model->tolerance, (int) REPORT, workdata, &fncount,
        &grcount, &fail);
  vmaxset(vmax);

  Memcpy(model->coef, coef, dm->p);

  Free(mask); Free(coef); Free(gradient); Free(workdata);
}

static double
nls_objective(int p, double *coef, void *workdata)
{
  SSIM_WKDAT wk = (SSIM_WKDAT) workdata;
  double val;
  int do_derivatives = 0;

  objective_SSIM(p, coef, &val, DNULLP, do_derivatives, wk);
  return val;
}

void
nls_gradient(int p, double *coef, double *gradient, void *workdata)
{
  SSIM_WKDAT wk = (SSIM_WKDAT) workdata;
  double val;
  int do_derivatives = 1;

  objective_SSIM(p, coef, &val, gradient, do_derivatives, wk);
  Memcpy(wk->nls_der, gradient, p);
}

void
objective_SSIM(int p, double *coef, double *fnc, double *gradient, int do_derivatives,
  SSIM_WKDAT wk)
{
  DIMS dm = wk->dm;
  TERMS terms = wk->terms;
  int job;
  double accum = 0.0, logdet = 0.0, factor, trace, varfnc, val, wts;
  double *score, *working, *u;

  job = 1;
  model_frame(terms, dm, coef, &(wk->phi), wk->x, wk->resid, wk->weights, job);

  varfnc = SQR(wk->phi) * (SQR(wk->phi) - 1.0);
  for (int i = 0; i < dm->n; i++) {
    logdet += log((wk->weights)[i]);
    accum  += SQR((wk->resid)[i]) / (wk->weights)[i];
  }
  val  = .5 * logdet + .5 * accum / varfnc;
  *fnc = val;

  if (!do_derivatives)
    return;

  /* get derivatives */
  score   = (double *) Calloc(dm->p, double);
  working = (double *) Calloc(dm->n, double);
  u       = (double *) Calloc(dm->n * dm->p, double);

  /* transformed model matrix and working residuals */
  for (int i = 0; i < dm->n; i++) {
    wts = sqrt((wk->weights)[i]);
    working[i] = (wk->resid)[i] / wts;
    for (int j = 0; j < dm->p; j++)
      u[i + j * dm->n] = (wk->x)[i + j * dm->n] / wts;
  }
  crossprod(score, u, dm->n, dm->n, dm->p, working, dm->n, dm->n, 1);
  factor = wk->phi / varfnc;
  scale_vec(-factor, score, 1, dm->p);

  /* 2nd term of the derivative with respect to coefficients */
  for (int j = 0; j < dm->p; j++) {
    accum = trace = 0.0;
    for (int i = 0; i < dm->n; i++) {
      wts = (wk->weights)[i];
      trace += (wk->x)[i + j * dm->n] / sqrt(wts);
      accum += SQR((wk->resid)[i] / wts) * (wk->x)[i + j * dm->n];
    }
    score[j] += trace - accum / varfnc;
  }

  Memcpy(gradient, score, p);

  Free(score); Free(working); Free(u);
}

void
update_phi(DIMS dm, double *phi, double *z, double *f, double tol)
{
  double conv, upper_phi;
  const double c = (1. + sqrt(5.)) * .5;
  PL_WKDAT workdata;

  workdata = (PL_WKDAT) Calloc(1, PL_workdata);

  /* info for a pseudo-likelihood object */
  workdata->dm  = dm;
  workdata->z   = z;
  workdata->f   = f;
  workdata->phi = *phi;

  /* call optimizer */
  upper_phi = *phi;
  do {
    *phi = brent(1., upper_phi, PL_objective, workdata, tol);
    conv = fabs(*phi - upper_phi);
    upper_phi *= c;
  } while (conv < tol);

  Free(workdata);
}

static double
PL_objective(double phi, void *workdata)
{ /* pseudo-likelihood objective function (to be invoked by brent's procedure)*/
  PL_WKDAT wk = (PL_WKDAT) workdata;
  DIMS dm = wk->dm;
  double accum = 0.0, resid, term, RSS;

  for (int i = 0; i < dm->n; i++) {
    resid  = (wk->z)[i] - phi * (wk->f)[i];
    term   = SQR((wk->f)[i]) * SQR(phi) * (SQR(phi) - 1.0);
    accum += SQR(resid) / term + log(term);
  }
  RSS = 0.5 * accum;
  wk->fnc = RSS;

  return RSS;
}

void
model_frame(TERMS terms, DIMS dm, double *coef, double *scale, double *x, double *resid,
  double *weights, int do_model_matrix)
{ /* compute residuals, weights and the local model matrix (if requested) */
  double alpha, beta, gamma, fnc;

  /* initialization */
  alpha = coef[0];
  beta  = coef[1];
  gamma = coef[2];

  /* residuals and weights */
  for (int i = 0; i < dm->n; i++) {
    fnc = R_pow((terms->luminance)[i], alpha) * R_pow((terms->contrast)[i], beta) *
          R_pow((terms->structure)[i], gamma);
    resid[i] = (terms->z)[i] - *scale * fnc;
    weights[i] = SQR(fnc);
  }

  if (!do_model_matrix)
    return;

  /* local model matrix */
  for (int i = 0; i < dm->n; i++) {
    x[i + 0 * dm->n] = alpha * R_pow((terms->luminance)[i], alpha - 1.) *
                       R_pow((terms->contrast)[i], beta) * R_pow((terms->structure)[i], gamma);
    x[i + 1 * dm->n] = beta * R_pow((terms->contrast)[i], beta - .1) *
                       R_pow((terms->luminance)[i], alpha) * R_pow((terms->structure)[i], gamma);
    x[i + 2 * dm->n] = gamma * R_pow((terms->structure)[i], gamma - 1.) *
                       R_pow((terms->luminance)[i], alpha) * R_pow((terms->contrast)[i], beta);
  }
}

void
mean_surface(TERMS terms, DIMS dm, double *coef, double *f)
{ /* evaluate the expected response function */
  double alpha, beta, gamma;

  /* initialization */
  alpha = coef[0];
  beta  = coef[1];
  gamma = coef[2];

  /* mean function (SSIM) */
  for (int i = 0; i < dm->n; i++) {
    f[i] = R_pow((terms->luminance)[i], alpha) * R_pow((terms->contrast)[i], beta) *
           R_pow((terms->structure)[i], gamma);
  }
}

double
SSIM_logLik(DIMS dm, double *resid, double *weights, double *scale)
{ /* evaluate the log-likelihood function for a heteroscedastic regression model */
  double accum = 0.0, varfnc, val;

  varfnc = SQR(*scale) * (SQR(*scale) - 1.0);
  for (int i = 0; i < dm->n; i++)
    accum += log(*weights++);

  accum += norm_sqr(resid, 1, dm->n) / varfnc;

  val  = (double) dm->n * (M_LN_SQRT_2PI + 0.5 * log(varfnc));
  val += 0.5 * accum;

  return -val;
}
