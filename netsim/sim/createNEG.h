#ifndef CREATENEG_H
#define CREATENEG_H

#include "r_aux.h"
#include "globali.h"
#include "sample.h"

MATRICEi *createNEG1(MATRICEi *ris, const MATRICEi *m);

SEXP createNEG(SEXP m);

#endif
