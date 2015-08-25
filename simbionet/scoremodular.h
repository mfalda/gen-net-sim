#ifndef SCORE_H
#define SCORE_H

#include "r_aux.h"
#include "globali.h"

VETTOREd *score1(VETTOREd *Sc, const VETTOREi *S, const VETTOREd *ST, const VETTOREd *Freq, int n, const VETTOREd *toll);

SEXP score(SEXP S, SEXP ST, SEXP Freq, SEXP n, SEXP toll);

#endif
