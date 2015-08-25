#ifndef CSF_H
#define CSF_H

#include "r_aux.h"
#include "globali.h"
#include "hist.h"
#include "sampleB.h"
#include "score_sf.h"
#include "scoremodular.h"

LISTA *connectivity_scalefree1(LISTA *ris, int N, int max_con, double gamma, double r_tol, double a_tol, double weight_mean, double weight_sd);

SEXP connectivity_scalefree(SEXP N, SEXP max_con, SEXP gamma, SEXP r_tol, SEXP a_tol, SEXP weight_mean, SEXP weight_sd);

#endif
