#ifndef TRPOL2_H
#define TRPOL2_H

#include <stdio.h>

#include "r_aux.h"
#include "globali.h"

VETTOREd *trpol2(VETTOREd *ris, int n, double x);

SEXP ret_trpol2(SEXP n, SEXP x);

#endif

