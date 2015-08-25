#ifndef ASSIGNNODES_UND_H
#define ASSIGNNODES_UND_H

#include "r_aux.h"
#include "sampleB.h"


VETTOREi *assign_nodes2_und(VETTOREi *ris, const MATRICEi *M, MATRICEi *Mdiscr, VETTOREi *h, const VETTOREi *hubs, const MATRICEd *Sc, const VETTOREi *Sin, int max_con);

SEXP assign_nodes_und(SEXP M, SEXP Mdiscr, SEXP h, SEXP hubs, SEXP Sc, SEXP Sin, SEXP max_con);

#endif
