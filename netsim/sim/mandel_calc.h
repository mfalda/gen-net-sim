#ifndef MANDELC_H
#define MANDELC_H

#include "r_aux.h"
#include "muParserDLL.h"

void OnError(muParserHandle_t hParser);

MATRICEi *mandel1(muParserHandle_t hParser_real, muParserHandle_t hParser_imag, MATRICEi *ris, VETTOREd *param, GString *reale, GString *img, double ris_x, double ris_y, int n);

SEXP mandel(SEXP param, SEXP real, SEXP img, SEXP n);

#endif

