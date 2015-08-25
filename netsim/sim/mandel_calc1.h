#ifndef MANDELC_H
#define MANDELC_H

#include "r_aux.h"
#include "pars.h"

MATRICEi *mandel1(MATRICEi *ris, VETTOREd *param, GString *reale, GString *img, double ris_x, double ris_y, int n);

SEXP mandel(SEXP param, SEXP real, SEXP img, SEXP n);

#endif

