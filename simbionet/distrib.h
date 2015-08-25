#ifndef DISTRIB_H
#define DISTRIB_H

#include "r_aux.h"

#define MAX_RIGA 1024


VETTOREd *leggi_seq_d(VETTOREd *ris, char *buf, int n);

VETTOREi *leggi_seq_i(VETTOREi *ris, char *buf, int n);

VETTOREd *rnorm_s(VETTOREd *ris, int n, double mean, double sd, const char *chi);

VETTOREd *rlnorm_s(VETTOREd *ris, int n, double mean, double sd, const char *chi);

VETTOREd *runif_s(VETTOREd *ris, int n, double a, double b, const char *chi);


#endif