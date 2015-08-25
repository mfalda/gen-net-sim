#ifndef SAMPLE_H
#define SAMPLE_H

// Questo file contiene algoritmi di R (sample.c)
#ifdef WIN32
#include <mem.h>
#endif

#include "r_aux.h"
#include "distrib.h"
#include "globali.h"

#ifdef MDEBUG
	#define sample_p(A, B, C, D, E, F) _sample_p(A, B, C, D, E, F, #A, __FILE__, __LINE__)
#else
	#define sample_p(A, B, C, D, E, F) _sample_p(A, B, C, D, E, F)
#endif
#ifdef MDEBUG
	VETTOREi *_sample_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p, const char *chi, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi *_sample_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p, const char *chi);
#endif

#ifdef MDEBUG
	#define sample(A, B, C, D, F) _sample(A, B, C, D, F, #A, __FILE__, __LINE__)
#else
	#define sample(A, B, C, D, F) _sample(A, B, C, D, F)
#endif
#ifdef MDEBUG
	VETTOREi *_sample(VETTOREi *ris, const VETTOREi *x, int k, int replace, const char *chi, const char *nome, const char *nomefile, int linea);
#else
	VETTOREi *_sample(VETTOREi *ris, const VETTOREi *x, int k, int replace, const char *chi);
#endif

#endif

