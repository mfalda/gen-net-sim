#ifndef CM_H
#define CM_H

#include "r_aux.h"
#include "globali.h"
#include "sampleB.h"
#include "module1.h"
#include "module2.h"
#include "module3.h"
#include "cluster_coeff.h"
#include "assign_nodes.h"

LISTA *connectivity_modular1(LISTA *ris, int N, int max_con, double gamma, GString *INdegree, double Cf_cl, VETTOREi *num_subnet, double r_tol, double a_tol, double weight_mean, double weight_sd);

SEXP connectivity_modular(SEXP N, SEXP max_con, SEXP gamma, SEXP INdegree, SEXP Cf_cl, SEXP num_subnet, SEXP r_tol, SEXP a_tol, SEXP weight_mean, SEXP weight_sd);

#endif
