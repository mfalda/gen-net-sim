#ifndef SCORE_SF_H
#define SCORE_SF_H

#include "r_aux.h"
#include "globali.h"

VETTOREd *score_sf1(VETTOREd *ris, const VETTOREi *S, const VETTOREd *ST, const VETTOREd *Freq, int n, const VETTOREd *toll);

SEXP score_sf(SEXP S, SEXP ST, SEXP Freq, SEXP n, SEXP toll);

#endif
