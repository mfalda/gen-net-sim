#ifndef SORT_H
#define SORT_H

#include "r_aux.h"

#define MIN_MERGESORT_LIST_SIZE    32

#ifdef MDEBUG
	#define ordine(A, B) _ordine(A, B, #A, __FILE__, __LINE__)
#else
	#define ordine(A, B) _ordine(A, B)
#endif
VETTOREi *ordine(VETTOREi *ris, VETTOREi *v);

#endif