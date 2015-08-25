#ifndef HMM_H
#define HMM_H

#include "r_aux.h"
#include "probmod.h"
#include "cluster_coeff.h"
#include "assign_nodes.h"
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

MATRICEi *HMM1(MATRICEi *ris, int N, double Cf_cl, double gamma, const VETTOREd *degree, const GString *INdegree, const Mod *modules, int L, const VETTOREd *prior_p_subnet1, int max_con, bool feedback, bool zero_nodes_with_indegree0, bool sepgraph, double r_tol, double a_tol, int iter);

SEXP HMM(SEXP N, SEXP Cf_cl, SEXP gamma, SEXP degree, SEXP INdegree, SEXP modules, SEXP prior_p_subnet, SEXP max_con, SEXP feedback, SEXP zero_nodes_with_indegree0, SEXP sepgraph, SEXP r_tol, SEXP a_tol, SEXP iter);

#endif
