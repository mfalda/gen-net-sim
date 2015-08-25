#ifndef MODULE2_H
#define MODULE2_H

#include "r_aux.h"
#include "globali.h"
#include "sampleB.h"
#include "mod2.h"
#include "probmod.h"
#include "triangola.h"

void module22(int num2, int Nrim, const MATRICEi *Mdiscr, const VETTOREi *h, const VETTOREi *h_new, const VETTOREi *Sin, const VETTOREi *Sout, const VETTOREd *STin, const VETTOREd *STout, const VETTOREd *Freq_in, const VETTOREd *Freq_out, int max_con, double Cf_c, const VETTOREd *toll, struct RisProbM *ris);

LISTA *module21(LISTA *ris, int num2, int Nrim, MATRICEi *Mdiscr, VETTOREi *h, VETTOREi *h_new, VETTOREi *Sin, VETTOREi *Sout, VETTOREd *STin, VETTOREd *STout, VETTOREd *Freq_in, VETTOREd *Freq_out, int max_con, double Cf_c, VETTOREd *toll);

SEXP module2(SEXP num1, SEXP Nrim, SEXP Mdiscr, SEXP h, SEXP h_new, SEXP Sin, SEXP Sout, SEXP STin, SEXP STout, SEXP Freq_in, SEXP Freq_out, SEXP max_con, SEXP Cf_c, SEXP toll);

#endif
