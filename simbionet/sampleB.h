#ifndef SAMPLEB_H
#define SAMPLEB_H

#include <stddef.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <Rmath.h>		/* for rxxx functions */
#include <errno.h>

#include "r_aux.h"
#include "sample.h"

// replace e' falso per default
#ifdef MDEBUG
	#define sampleB_p(A, B, C, D, E) _sampleB_p(A, B, C, D, E, #A, __FILE__, __LINE__)
#else
	#define sampleB_p(A, B, C, D, E) _sampleB_p(A, B, C, D, E)
#endif
#ifdef MDEBUG
	VETTOREi *_sampleB_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi *_sampleB_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p);
#endif

#ifdef MDEBUG
	#define sampleB(A, B, C, D) _sampleB(A, B, C, D, #A, __FILE__, __LINE__)
#else
	#define sampleB(A, B, C, D) _sampleB(A, B, C, D)
#endif
// replace e' falso per default
#ifdef MDEBUG
	VETTOREi *_sampleB(VETTOREi *ris, const VETTOREi *x, int k, int replace, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi *_sampleB(VETTOREi *ris, const VETTOREi *x, int k, int replace);
#endif


#endif

