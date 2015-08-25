#ifndef TARGET_H
#define TARGET_H

#include <math.h>

#include "r_aux.h"
#include "globali.h"
#include "ext_parser.h"
#include "boole_result.h"

//~ typedef double (FUNZIONE1) (double);

double calcola_f(int n, double t);

VETTOREd *target1(VETTOREd *ris, const VETTOREd *n, const LISTA *R, const MATRICEd *m, const MATRICEi *N, const MATRICEd *ext_in, const LISTA *ext_fun, double t, muParserHandle_t *hParsers, double sd_noise);

SEXP target(SEXP n, SEXP r, SEXP m, SEXP N, SEXP ext_in, SEXP ext_fun, SEXP t);

#endif
