#ifndef HIST_H
#define HIST_H

// Questo file contiene un algoritmo tratto da R (binning.c)

#include "r_aux.h"

#ifdef MDEBUG
	#define hist1(A, B, C, D, E, F) _hist1(A, B, C, D, E, F, #A, __FILE__, __LINE__)
#else
	#define hist1(A, B, C, D, E, F) _hist1(A, B, C, D, E, F)
#endif
#ifdef MDEBUG
	VETTOREi *_hist1(VETTOREi *ris, VETTOREi *x, VETTOREi *breaks, int right, int include_border, int naok, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi *_hist1(VETTOREi *ris, VETTOREi *x, VETTOREi *breaks, int right, int include_border, int naok);
#endif

SEXP hist(SEXP x, SEXP breaks, SEXP right, SEXP include_border, SEXP naok);

#endif
