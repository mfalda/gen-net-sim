#ifndef LSODA_H
#define LSODA_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "r_aux.h"

typedef struct {
	double m;
	double b;
	double k;
	double F;
	double o;
} Params;

MATRICEd *lsoda_oscillatore1(MATRICEd *ris, LISTA *parms, VETTOREd *X0, VETTOREd *times, const char *metodo, double atol, double rtol, double stat_thr, double stat_width);

SEXP lsoda_oscillatore(SEXP parms, SEXP X0, SEXP times, SEXP metodo, SEXP atol, SEXP rtol, SEXP stat_thr, SEXP stat_width);

#endif
