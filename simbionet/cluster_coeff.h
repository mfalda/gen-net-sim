#ifndef CLUSTCOEFF_H
#define CLUSTCOEFF_H

#include "r_aux.h"
#include "globali.h"

VETTOREd *cluster_coeff2(VETTOREd *Cg, const MATRICEi *W, double *coeff);

LISTA *cluster_coeff1(LISTA *ris, MATRICEi *W);

SEXP cluster_coeff(SEXP W);

#endif
