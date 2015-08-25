#ifndef ASSIGNNODES_H
#define ASSIGNNODES_H

#include "r_aux.h"
#include "globali.h"
#include "sampleB.h"


VETTOREi *assign_nodes2(VETTOREi *ris, const MATRICEi *M, MATRICEi *Mdiscr, VETTOREi *h, const VETTOREi *hubs, const MATRICEd *Sc, const VETTOREi *Sin, int max_con);

LISTA *assign_nodes1(LISTA * ris, MATRICEi *M, MATRICEi *Mdiscr, VETTOREi *h, VETTOREi *hubs, MATRICEd *Sc, VETTOREi *Sin, int max_con);

SEXP assign_nodes(SEXP M, SEXP Mdiscr, SEXP h, SEXP hubs, SEXP Sc, SEXP Sin, SEXP max_con);

#endif
