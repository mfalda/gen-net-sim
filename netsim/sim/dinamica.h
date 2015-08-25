#ifndef DINAMICA_H
#define DINAMICA_H

#include "r_aux.h"
#include "globali.h"
#include "target.h"

MATRICEd *dinamica1(MATRICEd *ris, LISTA *parms, LISTA *regole, LISTA *ext_fun, VETTOREd *n0, VETTOREd *times, muParserHandle_t *hParsers, double stat_thr, double stat_width, double sd_noise);

SEXP dinamica(SEXP parms, SEXP regole, SEXP ext_fun, SEXP n0, SEXP times);

#endif
