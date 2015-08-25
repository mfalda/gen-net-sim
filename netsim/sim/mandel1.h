#ifndef MANDEL1_H
#define MANDEL1_H

#include "r_aux.h"

#include <complex.h>

MATRICEi *mandel1(MATRICEi *ris, VETTOREd *param, double ris_x, double ris_y, int n);

SEXP mandel(SEXP param, SEXP n);

#endif

