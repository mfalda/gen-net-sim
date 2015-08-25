#ifndef SIMULATENET_H
#define SIMULATENET_H

#include "r_aux.h"
#include "globali.h"
#include "connectivity_random.h"
#include "connectivity_scalefree.h"
#include "connectivity_modular.h"
#include "connectivity_geometric.h"
#include "createRules.h"
#include "createNEG.h"
#include "write_table.h"
#include "dinamica.h"
#include "lsoda.h"


enum TParam { ESTERNO, UNIF, NORM, LOG_NORM };

enum NParam { LAMBDA, ALPHA, BETA, XMIN, XMAX, XZERO };

#define MAXPAR 6

typedef struct {
	LISTA *ris;
	LISTA *r;
} LISTA2;

LISTA2 *simulatenet1(LISTA2 *lst, int N, GString *connectivity, int max_reg, double gamma, GString *INdegree, double Cf_cl, VETTOREi *num_subnet, double kappa, GString *f_pr_and, VETTOREd *Xmin, VETTOREd *Xmax, VETTOREd *lambda, VETTOREd *x0, VETTOREd *weight_par, GString *act_fun, VETTOREd *alpha, VETTOREd *beta, VETTOREd *times, GString *method, int save, int ind_itera, int pad_reti, VETTOREd *ext1, VETTOREd *ext2, muParserHandle_t hParser, double stat_thr, double stat_width, int *stat);

SEXP simulatenet(SEXP N, SEXP connectivity, SEXP max_reg, SEXP gamma, SEXP INdegree, SEXP Cf_cl, SEXP num_subnet, SEXP kappa, SEXP f_pr_and, SEXP act_fun, SEXP alpha, SEXP beta, SEXP lambda, SEXP Xmin, SEXP Xmax, SEXP x0, SEXP weight_par, SEXP param, SEXP times, SEXP stat_thr, SEXP stat_width, SEXP method, SEXP num_reti, SEXP save);

#endif
