#ifndef INTERF_H
#define INTERF_H

#include "r_aux.h"
#include "globali.h"
#include "distrib.h" // per runif_s

VETTOREd *ret_vettore1(VETTOREd *ris, const VETTOREi *const vett, double num);
MATRICEd *ret_matrice1(MATRICEd *ris, const MATRICEd *matr);
LISTA *ret_lista1(LISTA *ris, const LISTA *lista);

SEXP ret_vettore(SEXP vett, SEXP num);
SEXP ret_matrice(SEXP matr);
SEXP ret_lista(SEXP lista);

#endif

