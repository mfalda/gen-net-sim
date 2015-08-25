#ifndef CLUSTCOEFF_H
#define CLUSTCOEFF_H

#include "r_aux.h"
#include "globali.h"

VETTOREd *cluster_coeff2(VETTOREd *ris, const MATRICEi *W, double *coeff);
LISTA *cluster_coeff1(LISTA *ris, const MATRICEi *W);

SEXP cluster_coeff(SEXP W);

#endif
