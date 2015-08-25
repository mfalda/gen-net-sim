#ifndef HMM_und_H
#define HMM_und_H

#include "r_aux.h"
#include "probmod_und.h"
#include "cluster_coeff.h"
#include "assign_nodes_und.h"
#include "check_conn.h"
#include "connetti_scalefree.h"
#include "sampleB.h"
#include "hist.h"
#include "scoremodular.h"
#include "write_table.h"

typedef struct {
	int codice;
	MATRICEi *rete;
	VETTOREi *hubs;
	int CC;
	bool autoreg;
	bool feedback;
	VETTOREi *hubsio;
	bool simm;
	int dim_m;
} Mod;

MATRICEi *HMM1_und(MATRICEi *ris, int N, double Cf_cl, double gamma, const VETTOREd *degree, const Mod *modules, int L, const VETTOREd *prior_p_subnet1, int max_con, bool sepgraph, double r_tol, double a_tol, int iter);

SEXP HMM_und(SEXP N, SEXP Cf_cl, SEXP gamma, SEXP degree, SEXP modules, SEXP prior_p_subnet, SEXP max_con, SEXP sepgraph, SEXP r_tol, SEXP a_tol, SEXP iter);

#endif
