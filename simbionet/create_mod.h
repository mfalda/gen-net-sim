#ifndef CREATEMOD_H
#define CREATEMOD_H

#include "r_aux.h"
#include "hubs.h"
#include "cluster_coeff.h"
#include "read_table.h"

typedef struct {
	int codice;
	MATRICEi *rete;
	VETTOREi *hubs;
	double CC;
	bool autoreg;
	bool feedback;
	VETTOREi *hubsio;
	bool simm;
	int dim_m;
} Mod;

GList *createMOD1(int m, bool auto1, int *len);

SEXP createMOD(SEXP m, SEXP auto1);

#endif
