#ifndef CONN_SF_H
#define CONN_SF_H

#include "r_aux.h"
#include "hist.h"
#include "scoremodular.h"
#include "sampleB.h"

MATRICEi *connetti_scalefree(MATRICEi *ris, const VETTOREd *STout, const VETTOREd *STin, const VETTOREd *dist, const VETTOREd *toll, int max_con, bool und/*=FALSE*/);

#endif
