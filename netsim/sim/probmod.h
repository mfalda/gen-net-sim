#ifndef PROBMOD_H
#define PROBMOD_H

#include <math.h>

#include "r_aux.h"
#include "globali.h"
#include "scoremodular.h"

struct RisProbM {
	double score;
	char label[3];
};

MATRICEd *probmod2(MATRICEd *ris, const MATRICEi *M, const VETTOREi *h, const VETTOREi *Sin, const VETTOREi *Sout, const VETTOREd *STin, const VETTOREd *STout, const VETTOREd *Freq_in, const VETTOREd *Freq_out, const VETTOREd *toll, struct RisProbM *ris1);

LISTA *probmod1(LISTA *ris, MATRICEi *M, VETTOREi *h, VETTOREi *Sin, VETTOREi *Sout, VETTOREd *STin, VETTOREd *STout, VETTOREd *Freq_in, VETTOREd *Freq_out, VETTOREd *toll);

SEXP probmod(SEXP M, SEXP h, SEXP Sin, SEXP Sout, SEXP STin, SEXP STout, SEXP Freq_in, SEXP Freq_out, SEXP toll);

#endif
