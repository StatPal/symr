/* $ID: brent.h, last updated 2018/07/23, F. Osorio */

#ifndef BRENT_H
#define BRENT_H

#include "base.h"

/* Brent's method for unidimensional optimization */
extern double brent(double, double, double (*f)(double, void *), void *, double);

#endif /* BRENT_H */
