#ifndef LSODA_H
#define LSODA_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "r_aux.h"
#include "globali.h"
#include "target.h"

typedef struct {
	int dim;
	VETTOREd *alpha;
	VETTOREd *theta;
	VETTOREd *k;
	LISTA *r;
	MATRICEd *M;
	MATRICEi *N;
	GString *af;
	MATRICEd *ei;
	LISTA *ef;
	VETTOREd *xmin;
	muParserHandle_t *hp;
	double sd;
} Params;

MATRICEd *lsoda1(MATRICEd *ris, LISTA *parms, LISTA *regole, LISTA *ext_fun, VETTOREd *X0, VETTOREd *times, const char *metodo, double atol, double rtol, muParserHandle_t *hParsers, double stat_thr, double stat_width, int *stat, double sd_noise);

SEXP lsoda(SEXP parms, SEXP regole, SEXP ext_fun, SEXP X0, SEXP times, SEXP metodo, SEXP atol, SEXP rtol);

#endif
