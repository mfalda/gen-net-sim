#ifndef LISTA_H
#define LISTA_H

#include <stdio.h>
#include <assert.h>
#include <stddef.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <Rmath.h>
#include <R_ext/Arith.h> // nan
#include <R_ext/Utils.h> // sort
#include <errno.h>

enum TIPO { INTERO, REALE, VETT, MATR, LST };

typedef struct  {
	double n_elem;
	enum TIPO *tipi;
	void **dati;
} LISTA;

// crea una lista con una matrice e due vettori
// list(m, v1, v2)
LISTA *lista3mvv_i(MATRICEi *m, VETTOREi *v1, VETTOREi *v2);

#endif
