#ifndef CR_H
#define CR_H

#include "r_aux.h"
#include "globali.h"
#include "sampleB.h"

int Fattoriale(int n);

LISTA *connectivity_random1(LISTA *ris, int N, int max_con, double k, double weight_mean, double weight_sd);

SEXP connectivity_random(SEXP N, SEXP max_con, SEXP k, SEXP weight_mean, SEXP weight_sd);

#endif
