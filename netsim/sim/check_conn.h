#ifndef CHECK_CONN_H
#define CHECK_CONN_H

#include "r_aux.h"
#include "globali.h"

VETTOREd *check_conn1(VETTOREd *ris, MATRICEi *Mdiscr);

SEXP check_conn(SEXP Mdiscr);

#endif
