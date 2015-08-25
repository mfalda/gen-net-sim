#ifndef CG_H
#define CG_H

#include "r_aux.h"
#include "globali.h"
#include "sampleB.h"

LISTA *connectivity_geometric1(LISTA *ris, int N, double k, double weight_mean, double weight_sd);

SEXP connectivity_geometric(SEXP N, SEXP k, SEXP weight_mean, SEXP weight_sd);

#endif
